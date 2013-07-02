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

(defn get-resource-list [edge-procs n-procs dims]                                                                                                                                           
  (let [procs (expt edge-procs dims)
         remainder (rem procs n-procs)
         rem-string (if (== 0 remainder) "" (str "+1:ppn=" remainder))]
    (str "-l nodes=" (quot procs n-procs) ":ppn=" n-procs rem-string)))

(defn submit-job [domain-size edge-procs executable results-directory q procs-per-node script-name dims time-steps wall-time]
  (clojure.java.shell/sh "qsub"
                         (str "-l walltime=0:" wall-time ":0")
                         (get-resource-list edge-procs procs-per-node dims)
                         (str "-vDOMAIN_SIZE=" domain-size ",PROCS_PER_EDGE=" edge-procs
                              ",PROG_NAME=" executable
                              ",RESULTS_DIR=" results-directory
                              ",NUM_STEPS=" time-steps)
                         (str "-N" executable domain-size "_" edge-procs)
                         (str "-q" q)
                         (str "-e" results-directory "/" executable "err_size" domain-size "procs" edge-procs)
                         script-name))

; ./submit-jobs.clj -s [40] -p [10 20 30 40 50] -e emsquare2d-one-line -n 12 -q janus-short -d /lustre/adhi1756/comm_only

(let [args-hash
      (first (cli *command-line-args*
           ["-s" "--domain-sizes" "Run the code with these domain sizes" :parse-fn #(read-string %)]
           ["-p" "--edge-procs" "Run the code with these many processors along the edges" :parse-fn #(read-string %)]
           ["-e" "--executable-name" "Use this executable"]
           ["-d" "--results-directory" "Write results to this folder"]
           ["-q" "--queue" "Submit to this queue"]
           ["-n" "--procs-per-node" "Number of processors per node on the cluster" :parse-fn #(read-string %)]
           ["-x" "--submission-script-name" "Name of the script to run"]
           ["-dim" "--dimensions" :parse-fn #(read-string %) :default 2]
           ["-t" "--time-steps" :parse-fn #(read-string %) :default [10000]]
           ["-w" "--wall-time" :parse-fn #(read-string %) :default [10]]))
      executable (if (contains? args-hash :executable-name) (:executable-name args-hash) "emsquare2d-strong")
      results-directory (if (contains? args-hash :results-directory)
                          (:results-directory args-hash)
                          ".")
      q (if (contains? args-hash :queue) (:queue args-hash) "lazy")
      procs-per-node (if (contains? args-hash :procs-per-node) (:procs-per-node args-hash) 8)
      time-steps (:time-steps args-hash)
      wall-time (:wall-time args-hash)
      script-name (if (contains? args-hash :submission-script-name) (:submission-script-name args-hash) "/projects/adhi1756/adi-prototype/qsemsquare2d-strong.q")
      ;; script-name "/projects/adhi1756/adi-prototype/qsemsquare2d-strong.q"
      dims (:dimensions args-hash)]
  (doseq [domain-size (:domain-sizes args-hash)
          edge-procs (:edge-procs args-hash)
          ts time-steps
          wt wall-time]
    (println (submit-job domain-size edge-procs executable results-directory q procs-per-node script-name dims ts wt))))