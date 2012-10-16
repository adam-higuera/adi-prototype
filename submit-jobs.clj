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

(defn get-resource-list [edge-procs]
  (if (== 0 (rem (expt edge-procs 2) 8))
    (str "-l nodes=" (quot (expt edge-procs 2) 7) ":ppn=7+1:ppn=" (rem (expt edge-procs 2) 7))
    (str "-l nodes=" (quot (expt edge-procs 2) 8) ":ppn=8+1:ppn=" (rem (expt edge-procs 2) 8))))

(defn submit-job [domain-size edge-procs executable results-directory q]
  (println (get-resource-list edge-procs))
  (clojure.java.shell/sh "qsub"
                         (get-resource-list edge-procs)
                         (str "-vDOMAIN_SIZE=" domain-size ",PROCS_PER_EDGE=" edge-procs
                              ",PROG_NAME=" executable
                              ",RESULTS_DIR=" results-directory)
                         (str "-N" executable domain-size "_" edge-procs)
                         (str "-q" q)
                         "/scr_verus/avh/Development/adi-prototype/qsemsquare2d-strong.q"))

(let [args-hash
      (cli *command-line-args*
           ["-s" "--domain-sizes" "Run the code with these domain sizes" :parse-fn #(read %)]
           ["-p" "--edge-procs" "Run the code with these many processors along the edges" :parse-fn #(read %)]
           ["-e" "--executable-name" "Use this executable"]
           ["-d" "--results-directory" "Write results to this folder"]
           ["-q" "--queue" "Submit to this queue"])
      executable (if (contains? args-hash :executable-name) (:executable-name args-hash) "emsquare2d-strong")
      results-directory (if (contains? args-hash :results-directory)
                          (:results-directory args-hash)
                          ".")
      q (if (contains? args-hash :queue) (:queue args-hash) "lazy")]
  (doseq [domain-size (:domain-sizes args-hash)
          edge-procs (:edge-procs args-hash)]
    (println (submit-job domain-size edge-procs executable results-directory))))
