#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
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
  yLine(world.split(world.rank() % procs_x + (world.rank() / (procs_x*procs_y)) * procs_x)),
  zLine(world.split(world.rank() % (procs_x*procs_y))),
  nSteps(n_steps),
  currentStep(0),
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
  procsX = xLine.size();
  procsY = yLine.size();
  procsZ = zLine.size();

  VacuumMatrixInitializer mat_init_x = VacuumMatrixInitializer(dx, dt, blockSize, determineBoundary(xLine));
  VacuumMatrixInitializer mat_init_y = VacuumMatrixInitializer(dy, dt, blockSize, determineBoundary(yLine));
  VacuumMatrixInitializer mat_init_z = VacuumMatrixInitializer(dz, dt, blockSize, determineBoundary(zLine));
  VacuumCouplingInitializer coupling_init_x = VacuumCouplingInitializer(& mat_init_x, blockSize, xLine);
  VacuumCouplingInitializer coupling_init_y = VacuumCouplingInitializer(& mat_init_y, blockSize, yLine);
  VacuumCouplingInitializer coupling_init_z = VacuumCouplingInitializer(& mat_init_z, blockSize, zLine);

  std::vector<AbstractMatrixInitializer*> mat_inits_x(blockSize*blockSize, & mat_init_x);
  std::vector<AbstractMatrixInitializer*> mat_inits_y(blockSize*blockSize, & mat_init_y);
  std::vector<AbstractMatrixInitializer*> mat_inits_z(blockSize*blockSize, & mat_init_z);
  std::vector<AbstractCouplingInitializer*> coupling_inits_x(blockSize*blockSize, & coupling_init_x);
  std::vector<AbstractCouplingInitializer*> coupling_inits_y(blockSize*blockSize, & coupling_init_y);
  std::vector<AbstractCouplingInitializer*> coupling_inits_z(blockSize*blockSize, & coupling_init_z);

  guardB = allocateGuardStorage();
  guardE = allocateGuardStorage();

  init->setOffsets(xLine, yLine, zLine);
  initFields(init);

  xUpdateRHSs = init->initCollection(mat_inits_x, coupling_inits_x, blockSize, xLine);
  yUpdateRHSs = init->initCollection(mat_inits_y, coupling_inits_y, blockSize, yLine);
  zUpdateRHSs = init->initCollection(mat_inits_z, coupling_inits_z, blockSize, zLine);

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

  for(unsigned int iy=0; iy < blockSize; iy++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      init->populateGuard(& guardE[0][0][(blockSize*iy + iz)*3],
			  & guardB[0][0][(blockSize*iy + iz)*3],
			  -1, iy, iz);
      init->populateGuard(& guardE[0][1][(blockSize*iy + iz)*3],
			  & guardB[0][1][(blockSize*iy + iz)*3],
			  blockSize, iy, iz);
    }
  }

  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      init->populateGuard(& guardE[1][0][(blockSize*ix + iz)*3],
			  & guardB[1][0][(blockSize*ix + iz)*3],
			  ix, -1, iz);
      init->populateGuard(& guardE[1][1][(blockSize*ix + iz)*3],
			  & guardB[1][1][(blockSize*ix + iz)*3],
			  ix, blockSize, iz);
    }
  }

  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy=0; iy < blockSize; iy++) {
      init->populateGuard(& guardE[2][0][(blockSize*ix + iy)*3],
			  & guardB[2][0][(blockSize*ix + iy)*3],
			  ix, iy, -1);
      init->populateGuard(& guardE[2][1][(blockSize*ix + iy)*3],
			  & guardB[2][1][(blockSize*ix + iy)*3],
			  ix, iy, blockSize);
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

  unsigned int steps_per_timing = 100;
  unsigned long* timings = new unsigned long[nSteps / steps_per_timing];

  std::cout << std::setprecision(10);
  for(currentStep=0; currentStep <= nSteps; currentStep++) {
    if (currentStep % dump_periodicity == 0 && n_dumps < total_dumps && dump) {
      std::ostringstream filename(std::ios::out);
      filename << dumpDir << "/dump3D_" << currentStep / dump_periodicity << "dt" << 10.0/nSteps << "dx" << dx << ".h5";
      this->dumpFields(filename.str());
      n_dumps++;
    }
    if (currentStep % steps_per_timing == 0 && world.rank() == 0) {
      gettimeofday(& t2, NULL);
      if(currentStep > 0)
	timings[currentStep / steps_per_timing - 1] =
	  1000000*(t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
	// std::cout << 1000000*(t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) << std::endl;
      t1=t2;
    }
    if(currentStep < nSteps)
      this->timeStep();
  }

  if(world.rank() == 0)
    dumpTimings(timings, nSteps / steps_per_timing, steps_per_timing);
}

void Simulation3D::dumpTimings(unsigned long* timings, hsize_t total_timings,
			       unsigned int steps_per_timing) {
  std::ostringstream filename(std::ios::out);
  filename << dumpDir << "/timing_s" << blockSize << "_p" << world.size() << ".h5";
  
  hid_t file_id=H5Fcreate(filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hid_t timingspace=H5Screate_simple(1, & total_timings, NULL);
  hid_t dx_space = H5Screate(H5S_SCALAR);
  hid_t dt_space = H5Screate(H5S_SCALAR);
  hid_t bs_space = H5Screate(H5S_SCALAR);
  hid_t ns_space = H5Screate(H5S_SCALAR);
  hid_t alg_name_space = H5Screate(H5S_SCALAR);
  hid_t atype = H5Tcopy(H5T_C_S1);
#ifndef YEE
  H5Tset_size(atype, xUpdateRHSs->getAlgName().length());
#else
  H5Tset_size(atype, std::string("Yee").length());
#endif
  H5Tset_strpad(atype, H5T_STR_NULLTERM);

  hid_t timing_dset_id = H5Dcreate(file_id, "timings", H5T_NATIVE_LONG, timingspace,
				   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t dx_attr_id = H5Acreate(timing_dset_id, "dx", H5T_NATIVE_DOUBLE, dx_space,
			       H5P_DEFAULT, H5P_DEFAULT);
  hid_t dt_attr_id = H5Acreate(timing_dset_id, "dt", H5T_NATIVE_DOUBLE, dt_space,
			       H5P_DEFAULT, H5P_DEFAULT);
  hid_t bs_attr_id = H5Acreate(timing_dset_id, "blockSize", H5T_NATIVE_UINT, bs_space,
			       H5P_DEFAULT, H5P_DEFAULT);
  hid_t ns_attr_id = H5Acreate(timing_dset_id, "stepsPerTiming", H5T_NATIVE_UINT, ns_space,
			       H5P_DEFAULT, H5P_DEFAULT);
  hid_t alg_attr_id = H5Acreate(timing_dset_id,"communicationStrategy", atype, alg_name_space,
				H5P_DEFAULT, H5P_DEFAULT);

  herr_t status = H5Dwrite(timing_dset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL,
			   H5P_DEFAULT, timings);
  status = H5Awrite(dx_attr_id, H5T_NATIVE_DOUBLE, & dx);
  status = H5Awrite(dt_attr_id, H5T_NATIVE_DOUBLE, & dt);
  status = H5Awrite(bs_attr_id, H5T_NATIVE_UINT, & blockSize);
  status = H5Awrite(ns_attr_id, H5T_NATIVE_UINT, & steps_per_timing);
#ifndef YEE
  status = H5Awrite(alg_attr_id, atype, xUpdateRHSs->getAlgName().c_str());
#else
  status = H5Awrite(alg_attr_id, atype, "Yee");
#endif

  H5Sclose(timingspace); H5Sclose(dx_space); H5Sclose(dt_space);
  H5Sclose(alg_name_space); H5Sclose(bs_space); H5Sclose(ns_space);
  H5Tclose(atype);
  H5Aclose(dx_attr_id); H5Aclose(dt_attr_id); H5Aclose(bs_attr_id);
  H5Aclose(alg_attr_id); H5Aclose(ns_attr_id);
  H5Dclose(timing_dset_id);
  H5Fclose(file_id);
}

void Simulation3D::dumpFields(std::string filename) {
  unsigned int p_x = xLine.rank(), p_y = yLine.rank(), p_z = zLine.rank();
  hsize_t start[4];
  start[0] = xLine.rank()*blockSize; start[1] = yLine.rank()*blockSize;
  start[2] = zLine.rank()*blockSize; start[3] = 0;

  hsize_t rhsx_start[3];
  rhsx_start[0] = p_y*blockSize; rhsx_start[1] = p_z*blockSize; rhsx_start[2] = p_x*blockSize;
  hsize_t rhsy_start[3];
  rhsy_start[0] = p_x*blockSize; rhsy_start[1] = p_z*blockSize; rhsy_start[2] = p_y*blockSize;
  hsize_t rhsz_start[3];
  rhsz_start[0] = p_x*blockSize; rhsz_start[1] = p_y*blockSize; rhsz_start[2] = p_z*blockSize;
  hsize_t rhs_count[3];
  rhs_count[0] = blockSize; rhs_count[1] = blockSize; rhs_count[2] = blockSize;

  hsize_t count[4];
  count[0] = blockSize; count[1] = blockSize; count[2] = blockSize; count[3] = 3;
  
  hsize_t dims[4];
  dims[0] = blockSize*procsX; dims[1] = blockSize*procsY; dims[2] = blockSize*procsZ; dims[3] = 3;
  hsize_t mem_dims[4];
  mem_dims[0] = blockSize; mem_dims[1] = blockSize; mem_dims[2] = blockSize; mem_dims[3] = 3;

  hid_t fa_p_list = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fa_p_list, world, MPI_INFO_NULL);

  hid_t file_id=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fa_p_list);
  H5Pclose(fa_p_list);

  hsize_t rhs_dims[3];
  rhs_dims[0] = blockSize*procsX; rhs_dims[1] = blockSize*procsY; rhs_dims[2] = blockSize*procsZ;
  hsize_t rhs_mem_dims[3];
  rhs_mem_dims[0] = blockSize; rhs_mem_dims[1] = blockSize; rhs_mem_dims[2] = blockSize;
  
  hid_t E_filespace = H5Screate_simple(4, dims, NULL);
  hid_t E_memspace = H5Screate_simple(4, mem_dims, NULL);
  hid_t B_filespace = H5Screate_simple(4, dims, NULL);
  hid_t B_memspace = H5Screate_simple(4, mem_dims, NULL);
  hid_t rhsx_filespace = H5Screate_simple(3, rhs_dims, NULL);
  hid_t rhsx_memspace = H5Screate_simple(3, rhs_mem_dims, NULL);
  hid_t rhsy_filespace = H5Screate_simple(3, rhs_dims, NULL);
  hid_t rhsy_memspace = H5Screate_simple(3, rhs_mem_dims, NULL);
  hid_t rhsz_filespace = H5Screate_simple(3, rhs_dims, NULL);
  hid_t rhsz_memspace = H5Screate_simple(3, rhs_mem_dims, NULL);

  hid_t E_dset_id = H5Dcreate(file_id, "E", H5T_NATIVE_DOUBLE, E_filespace,
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t E_wr_p_list = H5Pcreate(H5P_DATASET_XFER);
  // H5Pset_dxpl_mpio(E_wr_p_list, H5FD_MPIO_COLLECTIVE);
  H5Sselect_hyperslab(E_filespace, H5S_SELECT_SET, start, NULL, count, NULL);
  herr_t status = H5Dwrite(E_dset_id, H5T_NATIVE_DOUBLE, E_memspace, E_filespace,
			   E_wr_p_list, E);
  H5Dclose(E_dset_id);
  H5Sclose(E_filespace);
  H5Sclose(E_memspace);
  H5Pclose(E_wr_p_list);

  hid_t B_dset_id = H5Dcreate(file_id, "B", H5T_NATIVE_DOUBLE, B_filespace,
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t B_wr_p_list = H5Pcreate(H5P_DATASET_XFER);
  // H5Pset_dxpl_mpio(B_wr_p_list, H5FD_MPIO_COLLECTIVE);
  H5Sselect_hyperslab(B_filespace, H5S_SELECT_SET, start, NULL, count, NULL);
  status = H5Dwrite(B_dset_id, H5T_NATIVE_DOUBLE, B_memspace, B_filespace,
		    B_wr_p_list, B);
  H5Dclose(B_dset_id);
  H5Sclose(B_filespace);
  H5Sclose(B_memspace);
  H5Pclose(B_wr_p_list);

  hid_t rhsx_dset_id = H5Dcreate(file_id, "rhsx", H5T_NATIVE_DOUBLE, rhsx_filespace,
				 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t rhsx_wr_p_list = H5Pcreate(H5P_DATASET_XFER);
  // H5Pset_dxpl_mpio(rhsx_wr_p_list, H5FD_MPIO_COLLECTIVE);
  H5Sselect_hyperslab(rhsx_filespace, H5S_SELECT_SET, rhsx_start, NULL, rhs_count, NULL);
  status = H5Dwrite(rhsx_dset_id, H5T_NATIVE_DOUBLE, rhsx_memspace, rhsx_filespace,
		    rhsx_wr_p_list, rhsx);
  H5Dclose(rhsx_dset_id);
  H5Sclose(rhsx_filespace);
  H5Sclose(rhsx_memspace);
  H5Pclose(rhsx_wr_p_list);

  hid_t rhsy_dset_id = H5Dcreate(file_id, "rhsy", H5T_NATIVE_DOUBLE, rhsy_filespace,
				 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t rhsy_wr_p_list = H5Pcreate(H5P_DATASET_XFER);
  // H5Pset_dxpl_mpio(rhsy_wr_p_list, H5FD_MPIO_COLLECTIVE);
  H5Sselect_hyperslab(rhsy_filespace, H5S_SELECT_SET, rhsy_start, NULL, rhs_count, NULL);
  status = H5Dwrite(rhsy_dset_id, H5T_NATIVE_DOUBLE, rhsy_memspace, rhsy_filespace,
		    rhsy_wr_p_list, rhsy);
  H5Dclose(rhsy_dset_id);
  H5Sclose(rhsy_filespace);
  H5Sclose(rhsy_memspace);
  H5Pclose(rhsy_wr_p_list);

  hid_t rhsz_dset_id = H5Dcreate(file_id, "rhsz", H5T_NATIVE_DOUBLE, rhsz_filespace,
				 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t rhsz_wr_p_list = H5Pcreate(H5P_DATASET_XFER);
  // H5Pset_dxpl_mpio(rhsx_wr_p_list, H5FD_MPIO_COLLECTIVE);
  H5Sselect_hyperslab(rhsz_filespace, H5S_SELECT_SET, rhsz_start, NULL, rhs_count, NULL);
  status = H5Dwrite(rhsz_dset_id, H5T_NATIVE_DOUBLE, rhsz_memspace, rhsz_filespace,
		    rhsz_wr_p_list, rhsz);
  H5Dclose(rhsz_dset_id);
  H5Sclose(rhsz_filespace);
  H5Sclose(rhsz_memspace);
  H5Pclose(rhsz_wr_p_list);

  H5Fclose(file_id);
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
  std::ostringstream filename(std::ios::out);
#ifndef YEE
  implicitUpdateM();
  // filename << dumpDir << "/post_impM" << currentStep << ".h5";
  // dumpFields(filename.str());
  // filename.str("");
  explicitUpdateP();
  // filename << dumpDir << "/post_expP" << currentStep << ".h5";
  // dumpFields(filename.str());
  // filename.str("");
  implicitUpdateP();
  // filename << dumpDir << "/post_impP" << currentStep << ".h5";
  // dumpFields(filename.str());
  // filename.str("");
  explicitUpdateM();
  // filename << dumpDir << "/post_expM" << currentStep << ".h5";
  // dumpFields(filename.str());
#else
  yeeUpdate();
#endif
}

void Simulation3D::implicitUpdateM() {
  std::ostringstream filename(std::ios::out);

  populateRHSM();

  // filename << dumpDir << "/PreTD" << currentStep << ".h5";
  // dumpFields(filename.str());
  // filename.str("");

  xUpdateRHSs->doLines(rhs_ptrs_x);
  yUpdateRHSs->doLines(rhs_ptrs_y);
  zUpdateRHSs->doLines(rhs_ptrs_z);

  writeRHSM();

  // filename << dumpDir << "/PostTD" << currentStep << ".h5";
  // dumpFields(filename.str());
  // filename.str("");

  getGuardF(B, guardB, UPWARDS);
  implicitMSubstituteB();
  getGuardF(E, guardE, DOWNWARDS);

  // filename << dumpDir << "/post_Subst.h5";
  // dumpFields(filename.str());
}

void Simulation3D::implicitUpdateP() {
  std::ostringstream filename(std::ios::out);
  populateRHSP();

  // filename << dumpDir << "/PreTD" << currentStep << ".h5";
  // dumpFields(filename.str());
  // filename.str("");

  xUpdateRHSs->doLines(rhs_ptrs_x);
  yUpdateRHSs->doLines(rhs_ptrs_y);
  zUpdateRHSs->doLines(rhs_ptrs_z);

  writeRHSP();
  getGuardF(B, guardB, UPWARDS);
  implicitPSubstituteB();
  getGuardF(E, guardE, DOWNWARDS);
}

void Simulation3D::explicitUpdateP() {
  getGuardF(B, guardB, UPWARDS);
  
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
      }
    }
  }

  getGuardF(E, guardE, DOWNWARDS);
}

void Simulation3D::explicitUpdateM() {
  getGuardF(B, guardB, UPWARDS);
  
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
	if(iy + 1 < blockSize)
	  B[offset] -= preFactorY*(E[offset_yp1 + 2] - E[offset + 2]);
	else
	  B[offset] -= preFactorY*(guardEYH[(blockSize*ix + iz)*3 + 2] - E[offset + 2]);

	// B_y couples d E_x / dz
	if(iz + 1 < blockSize)
	  B[offset + 1] -= preFactorZ*(E[offset_zp1] - E[offset]);
	else
	  B[offset + 1] -= preFactorZ*(guardEZH[(blockSize*ix + iy)*3] - E[offset]);

	// B_z couples with d E_y / dx
	if(ix + 1 < blockSize)
	  B[offset + 2] -= preFactorX*(E[offset_xp1 + 1] - E[offset + 1]);
	else {
	  B[offset + 2] -= preFactorX*(guardEXH[(blockSize*iy + iz)*3 + 1] - E[offset + 1]);
	}
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
	  E[offset + 1] -= preFactorX*(tmp_field[offset + 2] - guardBXL[(blockSize*iy + iz)*3 + 2]);

	// d B_x / dy couples with E_z
	if(iy > 0)
	  E[offset + 2] -= preFactorY*(tmp_field[offset] - tmp_field[offset_ym1]);
	else
	  E[offset + 2] -= preFactorY*(tmp_field[offset] - guardBYL[(blockSize*ix + iz)*3]);

	// d B_y / dz couples with E_x
	if(iz > 0)
	  E[offset] -= preFactorZ*(tmp_field[offset + 1] - tmp_field[offset_zm1 + 1]);
	else
	  E[offset] -= preFactorZ*(tmp_field[offset + 1] - guardBZL[(blockSize*ix + iy)*3 + 1]);
      }
    }
  }

  getGuardF(E, guardE, DOWNWARDS);
}


void Simulation3D::implicitPSubstituteB() {
  double* guardBL = guardB[0][0]; double* guardBH = guardB[0][1];
  double* guardEH = guardE[0][1];
  for(unsigned int iy=0; iy < blockSize; iy++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int guard_offsetB = (blockSize*iy + iz)*3 + 1;
      unsigned int guard_offsetE = (blockSize*iy + iz)*3 + 2;
      // E_z coupled with d B_y / dx
      for(unsigned int ix=0; ix < blockSize; ix++) {
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3 + 1;
	unsigned int B_offset = (((ix-1)*blockSize + iy)*blockSize + iz)*3 + 1;
	unsigned int E_offset = (((ix)*blockSize + iy)*blockSize + iz)*3 + 2;
	if(ix > 0)
	  E[E_offset] += preFactorX*(B[B_offset2] - B[B_offset]);
	else
	  E[E_offset] += preFactorX*(B[B_offset2] - guardBL[guard_offsetB]);
      }
    }
  }

  guardBL = guardB[2][0]; guardBH = guardB[2][1]; guardEH = guardE[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy=0; iy < blockSize; iy++) {
      unsigned int guard_offsetB = (blockSize*ix + iy)*3;
      unsigned int guard_offsetE = (blockSize*ix + iy)*3 + 1;
      // E_y coupled with d B_x / dz
      for(unsigned int iz=0; iz < blockSize; iz++) {
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + (iz-1))*3;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 1;
	if(iz > 0)
	  E[E_offset] += preFactorZ*(B[B_offset2] - B[B_offset]);
	else
	  E[E_offset] += preFactorZ*(B[B_offset2] - guardBL[guard_offsetB]);
      }
    }
  }

  guardBL = guardB[1][0]; guardBH = guardB[1][1]; guardEH = guardE[1][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int guard_offsetB = (blockSize*ix + iz)*3 + 2;
      unsigned int guard_offsetE = (blockSize*ix + iz)*3;
      // E_x coupled with d B_z / dy
      for(unsigned int iy=0; iy < blockSize; iy++) {
	unsigned int B_offset = ((ix*blockSize + (iy-1))*blockSize + iz)*3 + 2;
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	if(iy > 0)
	  E[E_offset] += preFactorY*(B[B_offset2] - B[B_offset]);
	else
	  E[E_offset] += preFactorY*(B[B_offset2] - guardBL[guard_offsetB]);
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
      // E_y coupled with d B_z / dx
      for(unsigned int ix=0; ix < blockSize; ix++) {
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int B_offset = (((ix-1)*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int E_offset = (((ix)*blockSize + iy)*blockSize + iz)*3 + 1;
	if(ix > 0)
	  E[E_offset] -= preFactorX*(B[B_offset2] - B[B_offset]);
	else {
	  E[E_offset] -= preFactorX*(B[B_offset2] - guardBL[guard_offsetB]);
	}
      }
    }
  }

  guardBL = guardB[2][0]; guardBH = guardB[2][1]; guardEH = guardE[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy=0; iy < blockSize; iy++) {
      unsigned int guard_offsetB = (blockSize*ix + iy)*3 + 1;
      unsigned int guard_offsetE = (blockSize*ix + iy)*3;
      // E_x coupled with d B_y / dz
      for(unsigned int iz=0; iz < blockSize; iz++) {
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3 + 1;
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + (iz-1))*3 + 1;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	if(iz > 0)
	  E[E_offset] -= preFactorZ*(B[B_offset2] - B[B_offset]);
	else
	  E[E_offset] -= preFactorZ*(B[B_offset2] - guardBL[guard_offsetB]);
      }
    }
  }

  guardBL = guardB[1][0]; guardBH = guardB[1][1]; guardEH = guardE[1][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iz=0; iz < blockSize; iz++) {
      unsigned int guard_offsetB = (blockSize*ix + iz)*3;
      unsigned int guard_offsetE = (blockSize*ix + iz)*3 + 2;
      // E_z coupled with d B_x / dy
      for(unsigned int iy=0; iy < blockSize; iy++) {
	unsigned int B_offset = ((ix*blockSize + (iy-1))*blockSize + iz)*3;
	unsigned int B_offset2 = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	if(iy > 0)
	  E[E_offset] -= preFactorY*(B[B_offset2] - B[B_offset]);
	else
	  E[E_offset] -= preFactorY*(B[B_offset2] - guardBL[guard_offsetB]);
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
      // B_z coupled with d E_y / dx
      for(unsigned int ix=0; ix < blockSize; ix++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int E_offset = (((ix)*blockSize + iy)*blockSize + iz)*3 + 1;
	unsigned int E_offset2 = (((ix+1)*blockSize + iy)*blockSize + iz)*3 + 1;
	if(ix+1 < blockSize)
	  rhsx[rhs_offset + ix] = B[B_offset] - preFactorX*(E[E_offset2]-E[E_offset]);
	else {
	  rhsx[rhs_offset + ix] = B[B_offset] - preFactorX*(guardEH[guard_offset]-E[E_offset]);
	}
      }
    }
  }

  guardEH = guardE[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy=0; iy < blockSize; iy++) {
      unsigned int rhs_offset = blockSize*(blockSize*ix + iy);
      unsigned int guard_offset = (blockSize*ix + iy)*3;
      // B_y coupled with d E_x / dz
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
      // B_x coupled with d E_z / dy
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
      // B_y coupled with d E_z / dx
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
      // B_x coupled with d E_y / dz
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
      // B_z coupled with d E_x / dy
      for(unsigned int iy=0; iy < blockSize; iy++) {
	unsigned int B_offset = ((ix*blockSize + iy)*blockSize + iz)*3 + 2;
	unsigned int E_offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int E_offset2 = ((ix*blockSize + (iy+1))*blockSize + iz)*3;
	if(iy+1 < blockSize)
	  rhsy[rhs_offset + iy] = B[B_offset] + preFactorY*(E[E_offset2]-E[E_offset]);
	else {
	  rhsy[rhs_offset + iy] = B[B_offset] + preFactorY*(guardEH[guard_offset]-E[E_offset]);
	}
      }
    }
  }
}

void Simulation3D::writeRHSP() {
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

void Simulation3D::writeRHSM() {
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

// FIXME - THIS IS A VERY NAUGHTY FUNCTION.  YOU SHOULD BE ASHAMED
void Simulation3D::getGuardF(double* F, double*** guardF, commDir dir) {
  int p = xLine.rank();
  int n_p = xLine.size();

  unsigned int i_send = dir == UPWARDS ? blockSize - 1 : 0;
  int p_shift = dir == UPWARDS ? 1 : -1;
  double* guardDest = dir == UPWARDS ? guardF[0][0] : guardF[0][1];
  bool up = dir == UPWARDS;
  bool partner_up = ((p + p_shift >= 0) && !up) || ((p + p_shift < n_p) && up);
  bool partner_down = ((p - p_shift >= 0) && up) || ((p - p_shift < n_p) && !up);

  fillSendbufX(i_send, F);
  if(p % 2 == 0) {
    if(partner_up)
      xLine.send(p + p_shift, 0, guardSendbuf, 3*blockSize*blockSize);

    if(partner_down)
      xLine.recv(p - p_shift, 0, guardDest, 3*blockSize*blockSize);
    else if (F == B) {
      i_send = i_send == 0 ? blockSize - 1 : 0;
      fillSendbufX(i_send, F);
      std::copy(guardSendbuf, guardSendbuf + 3*blockSize*blockSize, guardDest);
    }
    else
      std::fill_n(guardDest, 3*blockSize*blockSize, 0);
  }
  else {
    if(partner_down)
      xLine.recv(p - p_shift, 0, guardDest, 3*blockSize*blockSize);
    else if (F == B) {
      i_send = i_send == 0 ? blockSize - 1 : 0;
      fillSendbufX(i_send, F);
      std::copy(guardSendbuf, guardSendbuf + 3*blockSize*blockSize, guardDest);
    }
    else
      std::fill_n(guardDest, 3*blockSize*blockSize, 0);

    if(partner_up)
      xLine.send(p + p_shift, 0, guardSendbuf, 3*blockSize*blockSize);
  }
 
  p = yLine.rank();
  n_p = yLine.size();
  i_send = dir == UPWARDS ? blockSize - 1 : 0;
  partner_up = ((p + p_shift >= 0) && !up) || ((p + p_shift < n_p) && up);
  partner_down = ((p - p_shift >= 0) && up) || ((p - p_shift < n_p) && !up);
  guardDest = dir == UPWARDS ? guardF[1][0] : guardF[1][1];

  fillSendbufY(i_send, F);
  if(p % 2 == 0) {
    if(partner_up)
      yLine.send(p + p_shift, 0, guardSendbuf, 3*blockSize*blockSize);

    if(partner_down)
      yLine.recv(p - p_shift, 0, guardDest, 3*blockSize*blockSize);
    else if (F == B) {
      i_send = i_send == 0 ? blockSize - 1 : 0;
      fillSendbufY(i_send, F);
      std::copy(guardSendbuf, guardSendbuf + 3*blockSize*blockSize, guardDest);
    }
    else
      std::fill_n(guardDest, 3*blockSize*blockSize, 0);
  }
  else {
    if(partner_down) {
      yLine.recv(p - p_shift, 0, guardDest, 3*blockSize*blockSize);
    }
    else if (F == B) {
      i_send = i_send == 0 ? blockSize - 1 : 0;
      fillSendbufY(i_send, F);
      std::copy(guardSendbuf, guardSendbuf + 3*blockSize*blockSize, guardDest);
    }
    else
      std::fill_n(guardDest, 3*blockSize*blockSize, 0);

    if(partner_up)
      yLine.send(p + p_shift, 0, guardSendbuf, 3*blockSize*blockSize);
  }

  p = zLine.rank();
  n_p = zLine.size();
  i_send = dir == UPWARDS ? blockSize - 1 : 0;
  partner_up = ((p + p_shift >= 0) && !up) || ((p + p_shift < n_p) && up);
  partner_down = ((p - p_shift >= 0) && up) || ((p - p_shift < n_p) && !up);
  guardDest = dir == UPWARDS ? guardF[2][0] : guardF[2][1];

  fillSendbufZ(i_send, F);
  if(p % 2 == 0) {
    if(partner_up)
      zLine.send(p + p_shift, 0, guardSendbuf, 3*blockSize*blockSize);

    if(partner_down)
      zLine.recv(p - p_shift, 0, guardDest, 3*blockSize*blockSize);
    else if (F == B) {
      i_send = i_send == 0 ? blockSize - 1 : 0;
      fillSendbufZ(i_send, F);
      std::copy(guardSendbuf, guardSendbuf + 3*blockSize*blockSize, guardDest);
    }
    else
      std::fill_n(guardDest, 3*blockSize*blockSize, 0);
  }
  else {
    if(partner_down)
      zLine.recv(p - p_shift, 0, guardDest, 3*blockSize*blockSize);
    else if (F == B) {
      i_send = i_send == 0 ? blockSize - 1 : 0;
      fillSendbufZ(i_send, F);
      std::copy(guardSendbuf, guardSendbuf + 3*blockSize*blockSize, guardDest);
    }
    else
      std::fill_n(guardDest, 3*blockSize*blockSize, 0);

    if(partner_up)
      zLine.send(p + p_shift, 0, guardSendbuf, 3*blockSize*blockSize);
  }
}

void Simulation3D::fillSendbufX(unsigned int ix, double* F) {
  for(unsigned int iy = 0; iy < blockSize; iy++) {   // Copy x-face into sendbuf
    for(unsigned int iz = 0; iz < blockSize; iz++) {
      for(unsigned int i = 0; i < 3; i++) {
	guardSendbuf[(blockSize*iy + iz)*3 + i] = F[(blockSize*(blockSize*ix + iy) + iz)*3 + i];
      }
    }
  }
}

void Simulation3D::fillSendbufY(unsigned int iy, double* F) {
  for(unsigned int ix = 0; ix < blockSize; ix++) {   // Copy y-face into sendbuf
    for(unsigned int iz = 0; iz < blockSize; iz++) {
      for(unsigned int i = 0; i < 3; i++)
	guardSendbuf[(blockSize*ix + iz)*3 + i] = F[(blockSize*(blockSize*ix + iy) + iz)*3 + i];
    }
  }
}

void Simulation3D::fillSendbufZ(unsigned int iz, double* F) {
  for(unsigned int ix = 0; ix < blockSize; ix++) {   // Copy z-face into sendbuf
    for(unsigned int iy = 0; iy < blockSize; iy++) {
      for(unsigned int i = 0; i < 3; i++)
	guardSendbuf[(blockSize*ix + iy)*3 + i] = F[(blockSize*(blockSize*ix + iy) + iz)*3 + i];
    }
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

void Simulation3D::yeeUpdate() {
  // update B
  double* guardEXH = guardE[0][1]; double* guardEYH = guardE[1][1]; double* guardEZH = guardE[2][1];
  for(unsigned int ix=0; ix < blockSize; ix++) {
    for(unsigned int iy = 0; iy < blockSize; iy++) {
      for(unsigned int iz = 0; iz < blockSize; iz++) {
	unsigned int offset = ((ix*blockSize + iy)*blockSize + iz)*3;
	unsigned int offset_xp1 = (((ix+1)*blockSize + iy)*blockSize + iz)*3;
	unsigned int offset_yp1 = ((ix*blockSize + (iy+1))*blockSize + iz)*3;
	unsigned int offset_zp1 = ((ix*blockSize + iy)*blockSize + (iz+1))*3;

	double E_x = E[offset]; double E_y = E[offset + 1]; double E_z = E[offset + 2];
	// B_x couples with dE_y / dz - d E_z / dy
	double E_z_yp1 = iy + 1 < blockSize ? E[offset_yp1 + 2] : guardEYH[(blockSize*ix + iz)*3 + 2];
	double E_y_zp1 = iz + 1 < blockSize ? E[offset_zp1 + 1] : guardEZH[(blockSize*ix + iy)*3 + 1];
	B[offset] += preFactorZ*(E_y_zp1 - E_y) - preFactorY*(E_z_yp1 - E_z);

	// B_y couples d E_z / dx - d E_x / dz
	double E_z_xp1 = ix + 1 < blockSize ? E[offset_xp1 + 2] : guardEXH[(blockSize*iy + iz)*3 + 2];
	double E_x_zp1 = iz + 1 < blockSize ? E[offset_zp1] : guardEZH[(blockSize*ix + iy)*3];
	B[offset + 1] += preFactorX*(E_z_xp1 - E_z) - preFactorZ*(E_x_zp1 - E_x);

	// B_z couples with d E_x / dy - d E_y / dx
	double E_x_yp1 = iy + 1 < blockSize ? E[offset_yp1] : guardEYH[(blockSize*ix + iz)*3];
	double E_y_xp1 = ix + 1 < blockSize ? E[offset_xp1 + 1] : guardEXH[(blockSize*iy + iz)*3 + 1];
	B[offset + 2] += preFactorY*(E_x_yp1 - E_x) - preFactorX*(E_y_xp1 - E_y);
      }
    }
  }
  // Share B
  getGuardF(B, guardB, UPWARDS);
  
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

	double B_x = B[offset]; double B_y = B[offset + 1]; double B_z = B[offset + 2];
	// d B_x / dz - d B_z / dx couples with E_y
	double B_z_xm1 = ix > 0 ? B[offset_xm1 + 2] : guardBXL[(blockSize*iy + iz)*3 + 2];
	double B_x_zm1 = iz > 0 ? B[offset_zm1] : guardBZL[(blockSize*ix + iy)*3];
	E[offset + 1] += preFactorZ*(B_x - B_x_zm1) - preFactorX*(B_z - B_z_xm1);

	// d B_y / d_x - d B_x / dy couples with E_z
	double B_x_ym1 = iy > 0 ? B[offset_ym1] : guardBYL[(blockSize*ix + iz)*3];
	double B_y_xm1 = ix > 0 ? B[offset_xm1 + 1] : guardBXL[(blockSize*iy + iz)*3 + 1];
	E[offset + 2] += preFactorX*(B_y - B_y_xm1) - preFactorY*(B_x - B_x_ym1);

	// d B_z / dy - d B_y / dz couples with E_x
	double B_y_zm1 = iz > 0 ? B[offset_zm1 + 1] : guardBZL[(blockSize*ix + iy)*3 + 1];
	double B_z_ym1 = iy > 0 ? B[offset_ym1 + 2] : guardBYL[(blockSize*ix + iz)*3 + 2];
	E[offset] += preFactorY*(B_z - B_z_ym1) - preFactorZ*(B_y - B_y_zm1);
      }
    }
  }
  // Share E
  getGuardF(E, guardE, DOWNWARDS);
}
