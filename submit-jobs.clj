#!/bin/bash lein-exec

(use '[leiningen.exec :only  (deps)])
(deps '[[org.clojure/math.numeric-tower "0.0.1"]
	[org.clojure/tools.cli "0.2.2"]])

(use 'clojure.math.numeric-tower)
(use '[clojure.tools.cli :only [cli]])
(require 'clojure.java.shell)

; (defn get-resource-list [edge-procs]
;   (if (<= edge-procs 8)
;     (str "-l nodes=" edge-procs ":ppn=" edge-procs)
;     (str "-l nodes=" (quot (expt edge-procs 2) 8) ":ppn=8+1:ppn=" (rem (expt edge-procs 2) 8))))

(defn get-resource-list [edge-procs n-procs]
  (if (== 0 (rem (expt edge-procs 2) n-procs))
    (str "-l nodes=" (quot (expt edge-procs 2) (- n-procs 1)) ":ppn=" (- n-procs 1)"+1:ppn=" (rem (expt edge-procs 2) (- n-procs 1)))
    (str "-l nodes=" (quot (expt edge-procs 2) n-procs) ":ppn=" n-procs "+1:ppn=" (rem (expt edge-procs 2) n-procs))))

(defn submit-job [domain-size edge-procs executable results-directory q procs-per-node]
  (clojure.java.shell/sh "qsub"
                         (get-resource-list edge-procs procs-per-node)
                         (str "-vDOMAIN_SIZE=" domain-size ",PROCS_PER_EDGE=" edge-procs
                              ",PROG_NAME=" executable
                              ",RESULTS_DIR=" results-directory)
                         (str "-N" executable domain-size "_" edge-procs)
                         (str "-q" q)
                         "/projects/adhi1756/adi-prototype/qsemsquare2d-strong.q"))

(let [args-hash
      (first (cli *command-line-args*
           ["-s" "--domain-sizes" "Run the code with these domain sizes" :parse-fn #(read-string %)]
           ["-p" "--edge-procs" "Run the code with these many processors along the edges" :parse-fn #(read-string %)]
           ["-e" "--executable-name" "Use this executable"]
           ["-d" "--results-directory" "Write results to this folder"]
           ["-q" "--queue" "Submit to this queue"]
           ["-n" "--procs-per-node" "Nubmer of processors per node on the cluster" :parse-fn #(read-string %)]))
      executable (if (contains? args-hash :executable-name) (:executable-name args-hash) "emsquare2d-strong")
      results-directory (if (contains? args-hash :results-directory)
                          (:results-directory args-hash)
                          ".")
      q (if (contains? args-hash :queue) (:queue args-hash) "lazy")
      procs-per-node (if (contains? args-hash :procs-per-node) (:procs-per-node args-hash) 8)]
  (doseq [domain-size (:domain-sizes args-hash)
          edge-procs (:edge-procs args-hash)]
    (println (submit-job domain-size edge-procs executable results-directory q procs-per-node))))
