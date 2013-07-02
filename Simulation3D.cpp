#include <fstream>
#include <sstream>
#include <iomanip>
#include "hdf5.h"
#include <boost/timer.hpp>
#include <sys/time.h>
#include "Simulation3D.hpp"
#define _USE_MATH_DEFINES

BoundaryLocation determineBoundary(mpi::communicator& world) {
  if(world.rank() == 0 && world.rank() == world.size()-1)
    return DOUBLE_BDY;
  else if(world.rank() == 0)
    return LOWER_BDY;
  else if(world.rank() == world.size() - 1)
    return UPPER_BDY;
  return NO_BDY;
}

Simulation3D::Simulation3D(double L_x, double L_y, double L_z,
			   double T,
			   unsigned int n_cells, unsigned int n_steps,
			   unsigned int procs_x, unsigned int procs_y, unsigned int procs_z,
			   unsigned int block_size,
			   std::string& dump_dir,
			   Simulation3DInitializer* init, mpi::communicator & world) :
  world(world),
  xLine(world.split(world.rank() / procs_x)),
  yLine(world.split(world.rank() % procs_x + world.rank() / (procs_x*procs_y))),
  zLine(world.split(world.rank() % (procs_x*procs_y))),
  nSteps(n_steps),
  dx(L_x/n_cells),
  dy(L_y/n_cells),
  dz(L_z/n_cells),
  dt(T/n_steps),
  blockSize(block_size),
  preFactorX(LIGHTSPEED*dt/(2*dx)),
  preFactorY(LIGHTSPEED*dt/(2*dy)),
  preFactorZ(LIGHTSPEED*dt/(2*dz)),
  E(new double[3*blockSize*blockSize*blockSize]),
  B(new double[3*blockSize*blockSize*blockSize]),
  tmp_field(new double[3*blockSize*blockSize*blockSize]),
  rhsx(new double[blockSize*blockSize*blockSize]),
  rhsy(new double[blockSize*blockSize*blockSize]),
  rhsz(new double[blockSize*blockSize*blockSize]),
  rhs_ptrs_x(new double*[blockSize*blockSize]),
  rhs_ptrs_y(new double*[blockSize*blockSize]),
  rhs_ptrs_z(new double*[blockSize*blockSize]),
  dumpDir(dump_dir)
{
  procsX = xLine.rank();
  procsY = yLine.rank();
  procsZ = zLine.rank();
  
  VacuumMatrixInitializer mat_init_x = VacuumMatrixInitializer(dx, dt, blockSize, determineBoundary(xLine));
  VacuumMatrixInitializer mat_init_y = VacuumMatrixInitializer(dy, dt, blockSize, determineBoundary(yLine));
  VacuumMatrixInitializer mat_init_z = VacuumMatrixInitializer(dz, dt, blockSize, determineBoundary(zLine));
  VacuumCouplingInitializer coupling_init_x = VacuumCouplingInitializer(& mat_init_x, blockSize, xLine);
  VacuumCouplingInitializer coupling_init_y = VacuumCouplingInitializer(& mat_init_y, blockSize, yLine);
  VacuumCouplingInitializer coupling_init_z = VacuumCouplingInitializer(& mat_init_z, blockSize, zLine);

  std::vector<AbstractMatrixInitializer*> mat_inits_x(blockSize, & mat_init_x);
  std::vector<AbstractMatrixInitializer*> mat_inits_y(blockSize, & mat_init_y);
  std::vector<AbstractMatrixInitializer*> mat_inits_z(blockSize, & mat_init_z);  
  std::vector<AbstractCouplingInitializer*> coupling_inits_x(blockSize, & coupling_init_x);
  std::vector<AbstractCouplingInitializer*> coupling_inits_y(blockSize, & coupling_init_y);
  std::vector<AbstractCouplingInitializer*> coupling_inits_z(blockSize, & coupling_init_z);

  init->setOffsets(xLine, yLine, zLine);
  initFields(init);
  
  xUpdateRHSs = init->initCollection(mat_inits_x, coupling_inits_x, blockSize, xLine);
  yUpdateRHSs = init->initCollection(mat_inits_y, coupling_inits_y, blockSize, yLine);
  zUpdateRHSs = init->initCollection(mat_inits_z, coupling_inits_z, blockSize, zLine);

  guardB = allocateGuardStorage();
  guardE = allocateGuardStorage();

  guardSendbuf = new double[3*blockSize*blockSize];
}

Simulation3D::~Simulation3D() {
  delete[] guardB[0][0];  delete[] guardB[0];   delete[] guardB;
  delete[] guardE[0][0];  delete[] guardE[0];   delete[] guardE;
  delete[] rhsx; delete[] rhsy; delete[] rhsz;
  delete[] tmp_field;
}

void Simulation3D::initFields(Simulation3DInitializer* init) {
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy=0; iy < blockSize; iy++) {
      for(unsigned int iz=0; iz < blockSize; iz++) {
	init->populate(E, B, ix, iy, iz);
      }
    }
  }

  for(unsigned int i=0; i < blockSize*blockSize; i++) {
    rhs_ptrs_x[i] = & rhsx[blockSize*i];
    rhs_ptrs_y[i] = & rhsy[blockSize*i];
    rhs_ptrs_z[i] = & rhsz[blockSize*i];
  }
}

void Simulation3D::simulate(bool dump, unsigned int dump_periodicity, unsigned int total_dumps) {
  unsigned int n_dumps = 0;
  timeval t1, t2;
  gettimeofday(& t1, NULL);

  std::cout << std::setprecision(10);
  for(unsigned int i=0; i <= nSteps; i++) {
    if (i % dump_periodicity == 0 && n_dumps < total_dumps && dump) {
      std::ostringstream filename(std::ios::out);
      filename << dumpDir << "/dump" << i / dump_periodicity << "dt" << 10.0/nSteps << "dx" << dx << ".h5";
      this->dumpFields(filename.str());
      n_dumps++;
    }
    if (i % 1000 == 0 && world.rank() == 0) {
      gettimeofday(& t2, NULL);
      if(i > 0)
	std::cout << 1000000*(t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) << std::endl;
      t1=t2;
    }
    this->timeStep();
  }
}

void Simulation3D::dumpFields(std::string filename) {
  hsize_t start[4];
  start[0] = xLine.rank()*blockSize; start[1] = yLine.rank()*blockSize;
  start[2] = zLine.rank()*blockSize; start[3] = 0;
  hsize_t count[4];
  count[0] = blockSize; count[1] = blockSize; count[2] = blockSize; count[3] = 3;
  
  hsize_t dims[4];
  dims[0] = blockSize*procsX; dims[1] = blockSize*procsY; dims[2] = blockSize*procsZ; dims[3] = 3;

  hid_t fa_p_list = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fa_p_list, world, MPI_INFO_NULL);

  hid_t file_id=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fa_p_list);
  H5Pclose(fa_p_list);

  hid_t E_filespace = H5Screate_simple(4, dims, NULL);
  hid_t B_filespace = H5Screate_simple(4, dims, NULL);
  hid_t E_dset_id = H5Dcreate(file_id, "E", H5T_NATIVE_DOUBLE, E_filespace,
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t B_dset_id = H5Dcreate(file_id, "B", H5T_NATIVE_DOUBLE, B_filespace,
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t E_wr_p_list = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(E_wr_p_list, H5FD_MPIO_COLLECTIVE);
  hid_t B_wr_p_list = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(B_wr_p_list, H5FD_MPIO_COLLECTIVE);

  H5Sselect_hyperslab(E_filespace, H5S_SELECT_SET, start, NULL, count, NULL);
  H5Sselect_hyperslab(B_filespace, H5S_SELECT_SET, start, NULL, count, NULL);

  herr_t status = H5Dwrite(E_dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
			   E_wr_p_list, E);
  status = H5Dwrite(B_dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    B_wr_p_list, B);

  H5Dclose(E_dset_id);
  H5Dclose(B_dset_id);
  H5Sclose(E_filespace);
  H5Sclose(B_filespace);
  H5Pclose(E_wr_p_list);
  H5Pclose(B_wr_p_list);
  H5Pclose(file_id);
}

double*** Simulation3D::allocateGuardStorage() {
  double* guardStorage = new double[6*3*blockSize*blockSize]; // 6 faces, 3 components, n^2
  double** innerIndexStorage = new double*[6]; // Pointers to storage for each face
  double*** guard = new double**[3]; // Pointers to the pair of pointers associated with each axis
  innerIndexStorage[0] = guardStorage;
  for (unsigned int i = 1; i < 6; i++) innerIndexStorage[i] = innerIndexStorage[i-1] + 3*blockSize*blockSize;
  for (unsigned int i = 0; i < 3; i++) guard[i] = & innerIndexStorage[2*i];
  return guard;
}

void Simulation3D::timeStep() {
  implicitUpdateM();
  explicitUpdateP();
  implicitUpdateP();
  explicitUpdateM();
}

void Simulation3D::implicitUpdateM() {
  populateRHSM();

  xUpdateRHSs->doLines(rhs_ptrs_x);
  yUpdateRHSs->doLines(rhs_ptrs_y);
  zUpdateRHSs->doLines(rhs_ptrs_z);

  writeRHSM();
  getGuardB();
  implicitMSubstituteB();
}

void Simulation3D::implicitUpdateP() {
  populateRHSP();

  xUpdateRHSs->doLines(rhs_ptrs_x);
  yUpdateRHSs->doLines(rhs_ptrs_y);
  zUpdateRHSs->doLines(rhs_ptrs_z);

  writeRHSP();
  getGuardB();
  implicitPSubstituteB();
}

void Simulation3D::explicitUpdateP() {
  getGuardB();
  
  // update B
  double* guardEXH = guardE[0][1]; double* guardEYH = guardE[1][1]; double* guardEZH = guardE[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy = 0; iy < blockSize; iy++) {
      for(unsigned int iz = 0; iz < blockSize; iz++) {
	unsigned int offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int offset_xp1 = (((ix+1)*blockSize + iy)*blockSize + iz)*3;
	unsigned int offset_yp1 = ((ix*blockSize + (iy+1))*blockSize + iz)*3;
	unsigned int offset_zp1 = ((ix*blockSize + iy)*blockSize + (iz+1))*3;

	for(unsigned int i=0; i < 3; i++) tmp_field[offset + i] = B[offset + i];

	// B_x couples with d E_y / dz
	if(iz + 1 < blockSize)
	  B[offset] += preFactorZ*(E[offset_zp1 + 1] - E[offset + 1]);
	else
	  B[offset] += preFactorZ*(guardEZH[(blockSize*ix + iy)*3 + 1] - E[offset + 1]);

	// B_y couples d E_z / dx
	if(ix + 1 < blockSize)
	  B[offset + 1] += preFactorX*(E[offset_xp1 + 2] - E[offset + 2]);
	else
	  B[offset + 1] += preFactorX*(guardEXH[(blockSize*iy + iz)*3 + 2] - E[offset + 2]);

	// B_z couples with d E_x / dy
	if(iy + 1 < blockSize)
	  B[offset + 2] += preFactorY*(E[offset_yp1] - E[offset]);
	else
	  B[offset + 2] += preFactorY*(guardEYH[(blockSize*ix + iz)*3] - E[offset]);
      }
    }
  }

  // update E
  double* guardBXL = guardB[0][0]; double* guardBYL = guardB[1][0]; double* guardBZL = guardB[2][0];
  double* guardBXH = guardB[0][1]; double* guardBYH = guardB[1][1]; double* guardBZH = guardB[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy = 0; iy < blockSize; iy++) {
      for(unsigned int iz = 0; iz < blockSize; iz++) {
	unsigned int offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int offset_xm1 = (((ix-1)*blockSize + iy)*blockSize + iz)*3;
	unsigned int offset_ym1 = ((ix*blockSize + (iy-1))*blockSize + iz)*3;
	unsigned int offset_zm1 = ((ix*blockSize + iy)*blockSize + (iz-1))*3;


	// d B_x / dz couples with E_y
	if(iz > 0)
	  E[offset + 1] += preFactorZ*(tmp_field[offset] - tmp_field[offset_zm1]);
	else
	  E[offset + 1] += preFactorZ*(tmp_field[offset] - guardBZL[(blockSize*ix + iy)*3]);

	// d B_y / dx couples with E_z
	if(ix > 0)
	  E[offset + 2] += preFactorX*(tmp_field[offset + 1] - tmp_field[offset_xm1 + 1]);
	else
	  E[offset + 2] += preFactorX*(tmp_field[offset + 1] - guardBXL[(blockSize*iy + iz)*3 + 1]);

	// d B_z / dy couples with E_x
	if(iy > 0)
	  E[offset] += preFactorY*(tmp_field[offset + 2] - tmp_field[offset_ym1 + 2]);
	else
	  E[offset] += preFactorY*(tmp_field[offset + 2] - guardBYL[(blockSize*ix + iz)*3 + 2]);

	if(ix + 1 == blockSize)
	  guardEXH[(blockSize*iy + iz)*3 + 2] += preFactorX*(guardBXH[(blockSize*iy + iz)*3 + 1] - tmp_field[offset_xm1 + 1]);
	if(iy + 1 == blockSize)
	  guardEYH[(blockSize*ix + iz)*3] += preFactorY*(guardBYH[(blockSize*ix + iz)*3 + 2] - tmp_field[offset_ym1 + 2]);
	if(iz + 1 == blockSize)
	  guardEZH[(blockSize*ix + iz)*3 + 1] += preFactorZ*(guardBZH[(blockSize*ix + iz)*3] - tmp_field[offset_zm1]);
      }
    }
  }
}

void Simulation3D::explicitUpdateM() {
  getGuardB();
  
  // update B
  double* guardEXH = guardE[0][1]; double* guardEYH = guardE[1][1]; double* guardEZH = guardE[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy = 0; iy < blockSize; iy++) {
      for(unsigned int iz = 0; iz < blockSize; iz++) {
	unsigned int offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int offset_xp1 = (((ix+1)*blockSize + iy)*blockSize + iz)*3;
	unsigned int offset_yp1 = ((ix*blockSize + (iy+1))*blockSize + iz)*3;
	unsigned int offset_zp1 = ((ix*blockSize + iy)*blockSize + (iz+1))*3;

	for(unsigned int i=0; i < 3; i++) tmp_field[offset + i] = B[offset + i];

	// B_x couples with d E_z / dy
	if(iz + 1 < blockSize)
	  B[offset] -= preFactorY*(E[offset_yp1 + 2] - E[offset + 2]);
	else
	  B[offset] -= preFactorY*(guardEYH[(blockSize*ix + iy)*3 + 2] - E[offset + 2]);

	// B_y couples d E_x / dz
	if(ix + 1 < blockSize)
	  B[offset + 1] -= preFactorZ*(E[offset_zp1] - E[offset]);
	else
	  B[offset + 1] -= preFactorZ*(guardEZH[(blockSize*iy + iz)*3] - E[offset]);

	// B_z couples with d E_y / dx
	if(iy + 1 < blockSize)
	  B[offset + 2] += preFactorX*(E[offset_xp1 + 1] - E[offset + 1]);
	else
	  B[offset + 2] += preFactorX*(guardEXH[(blockSize*ix + iz)*3 + 1] - E[offset + 1]);
      }
    }
  }

  // update E
  double* guardBXL = guardB[0][0]; double* guardBYL = guardB[1][0]; double* guardBZL = guardB[2][0];
  double* guardBXH = guardB[0][1]; double* guardBYH = guardB[1][1]; double* guardBZH = guardB[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy = 0; iy < blockSize; iy++) {
      for(unsigned int iz = 0; iz < blockSize; iz++) {
	unsigned int offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int offset_xm1 = (((ix-1)*blockSize + iy)*blockSize + iz)*3;
	unsigned int offset_ym1 = ((ix*blockSize + (iy-1))*blockSize + iz)*3;
	unsigned int offset_zm1 = ((ix*blockSize + iy)*blockSize + (iz-1))*3;


	// d B_z / dx couples with E_y
	if(ix > 0)
	  E[offset + 1] -= preFactorX*(tmp_field[offset + 2] - tmp_field[offset_xm1 + 2]);
	else
	  E[offset + 1] -= preFactorX*(tmp_field[offset] - guardBZL[(blockSize*ix + iy)*3]);

	// d B_x / dy couples with E_z
	if(iy > 0)
	  E[offset + 2] -= preFactorY*(tmp_field[offset] - tmp_field[offset_ym1]);
	else
	  E[offset + 2] -= preFactorY*(tmp_field[offset] - guardBXL[(blockSize*iy + iz)*3]);

	// d B_y / dz couples with E_x
	if(iz > 0)
	  E[offset] -= preFactorZ*(tmp_field[offset + 1] - tmp_field[offset_zm1 + 1]);
	else
	  E[offset] -= preFactorZ*(tmp_field[offset + 1] - guardBYL[(blockSize*ix + iz)*3 + 1]);

	if(ix + 1 == blockSize)
	  guardEXH[(blockSize*iy + iz)*3 + 1] -= preFactorX*(guardBXH[(blockSize*iy + iz)*3 + 2] - tmp_field[offset_xm1 + 2]);
	if(iy + 1 == blockSize)
	  guardEYH[(blockSize*ix + iz)*3 + 2] += preFactorY*(guardBYH[(blockSize*ix + iz)*3] - tmp_field[offset_ym1]);
	if(iz + 1 == blockSize)
	  guardEYH[(blockSize*ix + iz)*3] -= preFactorZ*(guardBZH[(blockSize*ix + iz)*3 + 1] - tmp_field[offset_zm1 + 1]);

      }
    }
  }
}


void Simulation3D::implicitPSubstituteB() {
  double* guardBL = guardB[0][0]; double* guardBH = guardB[0][1];
  double* guardEH = guardE[0][1];
  for(unsigned int iy=0; iy < blockSize; iy++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int guard_offsetB = (blockSize*iy + iz)*3 + 1;
      unsigned int guard_offsetE = (blockSize*iy + iz)*3 + 2;
      for(unsigned int ix=0; ix < blockSize; ix++) {
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3 + 1;
	unsigned int B_offset = (((ix-1)*blockSize + iy)*blockSize + iz)*3 + 1;
	unsigned int E_offset = (((ix)*blockSize + iy)*blockSize + iz)*3 + 2;
	if(ix > 0)
	  E[E_offset] += preFactorX*(B[B_offset2] - B[B_offset]);
	else
	  E[E_offset] += preFactorX*(B[B_offset2] - guardBL[guard_offsetB]);
	if(ix + 1 == blockSize)
	  guardEH[guard_offsetE] += preFactorX*(guardBH[guard_offsetB] - B[B_offset2]);
      }
    }
  }

  guardBL = guardB[2][0]; guardBH = guardB[2][1]; guardEH = guardE[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy=0; iy < blockSize; iy++) {
      unsigned int guard_offsetB = (blockSize*ix + iy)*3;
      unsigned int guard_offsetE = (blockSize*ix + iy)*3 + 1;
      for(unsigned int iz=0; iz < blockSize; iz++) {
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + (iz-1))*3;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 1;
	if(iz > 0)
	  E[E_offset] += preFactorZ*(B[B_offset2] - B[B_offset]);
	else
	  E[E_offset] += preFactorZ*(B[B_offset2] - guardBL[guard_offsetB]);
	if(ix + 1 == blockSize)
	  guardEH[guard_offsetE] += preFactorZ*(guardBH[guard_offsetB] - B[B_offset2]);
      }
    }
  }

  guardBL = guardB[1][0]; guardBH = guardB[1][1]; guardEH = guardE[1][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int guard_offsetB = (blockSize*ix + iz)*3 + 2;
      unsigned int guard_offsetE = (blockSize*ix + iz)*3;
      for(unsigned int iy=0; iy < blockSize; iy++) {
	unsigned int B_offset = ((ix*blockSize + (iy-1))*blockSize + iz)*3 + 2;
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int E_offset2 = ((ix*blockSize + (iy+1))*blockSize + iz)*3;
	if(iy > 0)
	  E[E_offset] += preFactorY*(B[B_offset2] - B[B_offset]);
	else
	  E[E_offset] += preFactorY*(B[B_offset2] - guardBL[guard_offsetB]);
	if(ix + 1 == blockSize)
	  guardEH[guard_offsetE] += preFactorY*(guardBH[guard_offsetB] - B[B_offset2]);
      }
    }
  }
}

void Simulation3D::implicitMSubstituteB() {
  double* guardBL = guardB[0][0]; double* guardBH = guardB[0][1];
  double* guardEH = guardE[0][1];
  for(unsigned int iy=0; iy < blockSize; iy++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int guard_offsetB = (blockSize*iy + iz)*3 + 2;
      unsigned int guard_offsetE = (blockSize*iy + iz)*3 + 1;
      for(unsigned int ix=0; ix < blockSize; ix++) {
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int B_offset = (((ix-1)*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int E_offset = (((ix)*blockSize + iy)*blockSize + iz)*3 + 1;
	if(ix > 0)
	  E[E_offset] -= preFactorX*(B[B_offset2] - B[B_offset]);
	else
	  E[E_offset] -= preFactorX*(B[B_offset2] - guardBL[guard_offsetB]);
	if(ix + 1 == blockSize)
	  guardEH[guard_offsetE] -= preFactorX*(guardBH[guard_offsetB] - B[B_offset2]);
      }
    }
  }

  guardBL = guardB[2][0]; guardBH = guardB[2][1]; guardEH = guardE[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy=0; iy < blockSize; iy++) {
      unsigned int guard_offsetB = (blockSize*ix + iy)*3 + 1;
      unsigned int guard_offsetE = (blockSize*ix + iy)*3;
      for(unsigned int iz=0; iz < blockSize; iz++) {
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3 + 1;
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + (iz-1))*3 + 1;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	if(iz > 0)
	  E[E_offset] -= preFactorZ*(B[B_offset2] - B[B_offset]);
	else
	  E[E_offset] -= preFactorZ*(B[B_offset2] - guardBL[guard_offsetB]);
	if(ix + 1 == blockSize)
	  guardEH[guard_offsetE] -= preFactorZ*(guardBH[guard_offsetB] - B[B_offset2]);
      }
    }
  }

  guardBL = guardB[1][0]; guardBH = guardB[1][1]; guardEH = guardE[1][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int guard_offsetB = (blockSize*ix + iz)*3 + 2;
      unsigned int guard_offsetE = (blockSize*ix + iz)*3;
      for(unsigned int iy=0; iy < blockSize; iy++) {
	unsigned int B_offset = ((ix*blockSize + (iy-1))*blockSize + iz)*3;
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	if(iy > 0)
	  E[E_offset] -= preFactorY*(B[B_offset2] - B[B_offset]);
	else
	  E[E_offset] -= preFactorY*(B[B_offset2] - guardBL[guard_offsetB]);
	if(ix + 1 == blockSize)
	  guardEH[guard_offsetE] -= preFactorY*(guardBH[guard_offsetB] - B[B_offset2]);
      }
    }
  }
}


void Simulation3D::populateRHSM() {
  double* guardEH = guardE[0][1];
  for(unsigned int iy=0; iy < blockSize; iy++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int rhs_offset = blockSize*blockSize*iy + blockSize*iz;
      unsigned int guard_offset = (blockSize*iy + iz)*3 + 1;
      for(unsigned int ix=0; ix < blockSize; ix++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int E_offset = (((ix)*blockSize + iy)*blockSize + iz)*3 + 1;
	unsigned int E_offset2 = (((ix+1)*blockSize + iy)*blockSize + iz)*3 + 1;
	if(ix+1 < blockSize)
	  rhsx[rhs_offset + ix] = B[B_offset] - preFactorX*(E[E_offset2]-E[E_offset]);
	else
	  rhsx[rhs_offset + ix] = B[B_offset] - preFactorX*(guardEH[guard_offset]-E[E_offset]);
      }
    }
  }

  guardEH = guardE[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy=0; iy < blockSize; iy++) {
      unsigned int rhs_offset = blockSize*(blockSize*ix + iy);
      unsigned int guard_offset = (blockSize*ix + iy)*3;
      for(unsigned int iz=0; iz < blockSize; iz++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 1;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int E_offset2 = ((ix*blockSize + iy)*blockSize + (iz+1))*3;
	if(iz+1 < blockSize)
	  rhsz[rhs_offset + iz] = B[B_offset] - preFactorZ*(E[E_offset2]-E[E_offset]);
	else
	  rhsz[rhs_offset + iz] = B[B_offset] - preFactorZ*(guardEH[guard_offset] -E[E_offset]);
      }
    }
  }

  guardEH = guardE[1][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int rhs_offset = blockSize*(blockSize*ix + iz);
      unsigned int guard_offset = (blockSize*ix + iz)*3 + 2;
      for(unsigned int iy=0; iy < blockSize; iy++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int E_offset2 = ((ix*blockSize + (iy+1))*blockSize + iz)*3 + 2;
	if(iy+1 < blockSize)
	  rhsy[rhs_offset + iy] = B[B_offset] - preFactorY*(E[E_offset2]-E[E_offset]);
	else
	  rhsy[rhs_offset + iy] = B[B_offset] - preFactorY*(guardEH[guard_offset]-E[E_offset]);
      }
    }
  }
}

void Simulation3D::populateRHSP() {
  double* guardEH = guardE[0][1];
  for(unsigned int iy=0; iy < blockSize; iy++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int rhs_offset = blockSize*blockSize*iy + blockSize*iz;
      unsigned int guard_offset = (blockSize*iy + iz)*3 + 2;
      for(unsigned int ix=0; ix < blockSize; ix++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 1;
	unsigned int E_offset = (((ix)*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int E_offset2 = (((ix+1)*blockSize + iy)*blockSize + iz)*3 + 2;
	if(ix+1 < blockSize)
	  rhsx[rhs_offset + ix] = B[B_offset] + preFactorX*(E[E_offset2]-E[E_offset]);
	else
	  rhsx[rhs_offset + ix] = B[B_offset] + preFactorX*(guardEH[guard_offset]-E[E_offset]);
      }
    }
  }

  guardEH = guardE[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy=0; iy < blockSize; iy++) {
      unsigned int rhs_offset = blockSize*(blockSize*ix + iy);
      unsigned int guard_offset = (blockSize*ix + iy)*3 + 1;
      for(unsigned int iz=0; iz < blockSize; iz++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 1;
	unsigned int E_offset2 = ((ix*blockSize + iy)*blockSize + (iz+1))*3 + 1;
	if(iz+1 < blockSize)
	  rhsz[rhs_offset + iz] = B[B_offset] + preFactorZ*(E[E_offset2]-E[E_offset]);
	else
	  rhsz[rhs_offset + iz] = B[B_offset] + preFactorZ*(guardEH[guard_offset] -E[E_offset]);
      }
    }
  }

  guardEH = guardE[1][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int rhs_offset = blockSize*(blockSize*ix + iz);
      unsigned int guard_offset = (blockSize*ix + iz)*3;
      for(unsigned int iy=0; iy < blockSize; iy++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int E_offset2 = ((ix*blockSize + (iy+1))*blockSize + iz)*3;
	if(iy+1 < blockSize)
	  rhsy[rhs_offset + iy] = B[B_offset] + preFactorY*(E[E_offset2]-E[E_offset]);
	else
	  rhsy[rhs_offset + iy] = B[B_offset] + preFactorY*(guardEH[guard_offset]-E[E_offset]);	  
      }
    }
  }
}

void Simulation3D::writeRHSP() { // FIXME !!!! once things compile
  for(unsigned int iy=0; iy < blockSize; iy++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int rhs_offset = blockSize*blockSize*iy + blockSize*iz;
      for(unsigned int ix=0; ix < blockSize; ix++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 1;
	B[B_offset] = rhsx[rhs_offset + ix];
      }
    }
  }

  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy=0; iy < blockSize; iy++) {
      unsigned int rhs_offset = blockSize*(blockSize*ix + iy);
      for(unsigned int iz=0; iz < blockSize; iz++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	B[B_offset] = rhsz[rhs_offset + iz];
      }
    }
  }

  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int rhs_offset = blockSize*(blockSize*ix + iz);
      for(unsigned int iy=0; iy < blockSize; iy++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	B[B_offset] = rhsy[rhs_offset + iy];
      }
    }
  }
}

void Simulation3D::writeRHSM() { // FIXME !!!! once things compile
  for(unsigned int iy=0; iy < blockSize; iy++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int rhs_offset = blockSize*blockSize*iy + blockSize*iz;
      for(unsigned int ix=0; ix < blockSize; ix++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	B[B_offset] = rhsx[rhs_offset + ix];
      }
    }
  }

  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy=0; iy < blockSize; iy++) {
      unsigned int rhs_offset = blockSize*(blockSize*ix + iy);
      for(unsigned int iz=0; iz < blockSize; iz++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 1;
	B[B_offset] = rhsz[rhs_offset + iz];
      }
    }
  }

  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int rhs_offset = blockSize*(blockSize*ix + iz);
      for(unsigned int iy=0; iy < blockSize; iy++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	B[B_offset] = rhsy[rhs_offset + iy];
      }
    }
  }
}


void Simulation3D::getGuardB() {
  unsigned int p = xLine.rank();
  unsigned int n_p = xLine.size();

  // even ranks exchange upwards
  unsigned int i_send = p % 2 == 0 ? blockSize - 1 : 0;
  double* guardDest = p % 2 == 0 ? guardB[0][0] : guardB[0][1];

  fillSendbufX(i_send);
  exchangeData(guardDest, guardSendbuf, xLine, UPWARDS);
  exchangeData(guardDest, guardSendbuf, xLine, DOWNWARDS);

  // even ranks exchange downwards
  i_send = p % 2 == 0 ? 0 : blockSize - 1;
  guardDest = p % 2 == 0 ? guardB[0][1] : guardB[0][0];

  fillSendbufX(i_send);
  exchangeData(guardDest, guardSendbuf, xLine, UPWARDS);
  exchangeData(guardDest, guardSendbuf, xLine, DOWNWARDS);

  p = yLine.rank();
  n_p = yLine.size();

  // even ranks exchange upwards
  i_send = p % 2 == 0 ? blockSize - 1 : 0;
  guardDest = p % 2 == 0 ? guardB[0][0] : guardB[0][1];

  fillSendbufY(i_send);
  exchangeData(guardDest, guardSendbuf, xLine, UPWARDS);
  exchangeData(guardDest, guardSendbuf, xLine, DOWNWARDS);

  // even ranks exchange downwards
  i_send = p % 2 == 0 ? 0 : blockSize - 1;
  guardDest = p % 2 == 0 ? guardB[0][1] : guardB[0][0];

  fillSendbufY(i_send);
  exchangeData(guardDest, guardSendbuf, xLine, UPWARDS);
  exchangeData(guardDest, guardSendbuf, xLine, DOWNWARDS);

  p = zLine.rank();
  n_p = zLine.size();

  // even ranks exchange upwards
  i_send = p % 2 == 0 ? blockSize - 1 : 0;
  guardDest = p % 2 == 0 ? guardB[0][0] : guardB[0][1];

  fillSendbufZ(i_send);
  exchangeData(guardDest, guardSendbuf, xLine, UPWARDS);
  exchangeData(guardDest, guardSendbuf, xLine, DOWNWARDS);

  // even ranks exchange downwards
  i_send = p % 2 == 0 ? 0 : blockSize - 1;
  guardDest = p % 2 == 0 ? guardB[0][1] : guardB[0][0];

  fillSendbufZ(i_send);
  exchangeData(guardDest, guardSendbuf, xLine, UPWARDS);
  exchangeData(guardDest, guardSendbuf, xLine, DOWNWARDS);  
}

void Simulation3D::fillSendbufX(unsigned int ix) {
  for(unsigned int iy = 0; iy < blockSize; iy++) {   // Copy x-face into sendbuf
    for(unsigned int iz = 0; iz < blockSize; iz++) {
      for(unsigned int i = 0; i < 3; i++)
	guardSendbuf[(blockSize*iy + iz)*3 + i] = B[(blockSize*(blockSize*ix + iy) + iz)*3 + i];
    }
  }
}

void Simulation3D::fillSendbufY(unsigned int iy) {
  for(unsigned int ix = 0; iy < blockSize; ix++) {   // Copy y-face into sendbuf
    for(unsigned int iz = 0; iz < blockSize; iz++) {
      for(unsigned int i = 0; i < 3; i++)
	guardSendbuf[(blockSize*ix + iz)*3 + i] = B[(blockSize*(blockSize*ix + iy) + iz)*3 + i];
    }
  }
}

void Simulation3D::fillSendbufZ(unsigned int iz) {
  for(unsigned int ix = 0; ix < blockSize; ix++) {   // Copy z-face into sendbuf
    for(unsigned int iy = 0; iy < blockSize; iy++) {
      for(unsigned int i = 0; i < 3; i++)
	guardSendbuf[(blockSize*ix + iz)*3 + i] = B[(blockSize*(blockSize*ix + iy) + iz)*3 + i];
    }
  }
}

void Simulation3D::exchangeData(double* guardStorage, double* sendbuf,
				mpi::communicator& world, commDir d) {
  switch(d) {
  case UPWARDS:
    if (world.rank() + 1 == world.size()) { // if no partner, symmetry BCs on B
      std::copy(guardSendbuf, guardSendbuf + 3*blockSize*blockSize, guardStorage);
      return;
    }
    world.send(world.rank() + 1, 0, sendbuf, 3*blockSize*blockSize);
    world.recv(world.rank() + 1, 0, guardStorage, 3*blockSize*blockSize);
    break;
  case DOWNWARDS:
    if (world.rank() == 0) { // if no partner, symmetry BCs on B
      std::copy(guardSendbuf, guardSendbuf + 3*blockSize*blockSize, guardStorage);
      return;
    }
    world.recv(world.rank() - 1, 0, guardStorage, 3*blockSize*blockSize);
    world.send(world.rank() - 1, 0, sendbuf, 3*blockSize*blockSize);
    break;
  }
}

Simulation3DInitializer::Simulation3DInitializer(double dx, double dy, double dz,
					     double L_x, double L_y, double L_z,
					     unsigned int blockSize) 
  : dx(dx), dy(dy), dz(dz),
    L_x(L_x), L_y(L_y), L_z(L_z),
    blockSize(blockSize) {}

void Simulation3DInitializer::setOffsets(const mpi::communicator & xLine,
					  const mpi::communicator & yLine,
					  const mpi::communicator & zLine) {
  x_offset=(dx*blockSize*xLine.rank());
  y_offset=(dy*blockSize*yLine.rank());
  z_offset=(dz*blockSize*zLine.rank());
}
