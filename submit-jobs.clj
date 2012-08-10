#!/bin/bash lein-exec

(use '[leiningen.exec :only  (deps)])
(deps '[[org.clojure/math.numeric-tower "0.0.1"]])

(use 'clojure.math.numeric-tower)
(require 'clojure.java.shell)

; (defn get-resource-list [edge-procs]
;   (if (<= edge-procs 8)
;     (str "-l nodes=" edge-procs ":ppn=" edge-procs)
;     (str "-l nodes=" (quot (expt edge-procs 2) 8) ":ppn=8+1:ppn=" (rem (expt edge-procs 2) 8))))

(defn get-resource-list [edge-procs]
  (if (== 0 (rem (expt edge-procs 2) 8))
    (str "-l nodes=" (quot (expt edge-procs 2) 7) ":ppn=7+1:ppn=" (rem (expt edge-procs 2) 7))
    (str "-l nodes=" (quot (expt edge-procs 2) 8) ":ppn=8+1:ppn=" (rem (expt edge-procs 2) 8))))

(defn submit-job [domain-size edge-procs]
  (println (get-resource-list edge-procs))
  (clojure.java.shell/sh "qsub"
                         (get-resource-list edge-procs)
                         (str "-vDOMAIN_SIZE=" domain-size ",PROCS_PER_EDGE=" edge-procs)
                         (str "-Nemsquare2d" domain-size "_" edge-procs)
                         "/scr_verus/avh/Development/adi-prototype/qsemsquare2d-strong.q"))

(doseq [domain-size [250 500 1000 2000]
        edge-procs [3 4 5 6 7 8 9 10]]
  (println (submit-job domain-size edge-procs)))

