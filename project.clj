(defproject kp "0.1.0-SNAPSHOT"
  :description "Kronig–Penney model CLI in Clojure"
  :url "https://example.com/kp"
  :license {:name "MIT"}
  :dependencies [[org.clojure/clojure "1.11.1"]
                 [org.clojure/tools.cli "1.1.230"]
                 [org.clojure/data.csv "1.0.1"]]
  :source-paths ["src"]
  :main kp.core
  :profiles {:uberjar {:aot :all}})
