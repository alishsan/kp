(ns kp.core-test
  (:require [clojure.test :refer :all]
            [kp.core :refer :all]))

(deftest test-parse-layer-token
  (testing "parse-layer-token function parses layer specifications"
    (is (= {:w 0.3 :V 0.0} (parse-layer-token "b:0.3")))
    (is (= {:w 0.2 :V 8.0} (parse-layer-token "U:8:0.2")))
    (is (= {:w 0.5 :V 0.0} (parse-layer-token "0.5")))))

(deftest test-parse-layers
  (testing "parse-layers function parses layer strings"
    (let [layers (parse-layers "b:0.3,U:8:0.2,b:0.3")]
      (is (= 3 (count layers)))
      (is (= {:w 0.3 :V 0.0} (first layers)))
      (is (= {:w 0.2 :V 8.0} (second layers)))
      (is (= {:w 0.3 :V 0.0} (nth layers 2))))))

(deftest test-generate-energy-grid
  (testing "generate-energy-grid function creates correct energy range"
    (let [energies (generate-energy-grid 0.0 10.0 5)]
      (is (= 6 (count energies)))
      (is (= 0.0 (first energies)))
      (is (= 10.0 (last energies)))
      (is (= 2.0 (nth energies 1)))
      (is (= 8.0 (nth energies 4))))))

(deftest test-cli-options
  (testing "cli-options contains expected options"
    (let [option-names (map first cli-options)]
      (is (some #{"-a"} option-names))
      (is (some #{"-b"} option-names))
      (is (some #{"-V"} option-names))
      (is (some #{"--layers"} option-names))
      (is (some #{"-o"} option-names)))))

(deftest test-usage
  (testing "usage function returns help string"
    (let [help-text (usage "test-summary")]
      (is (string? help-text))
      (is (.contains help-text "Kronig–Penney model CLI"))
      (is (.contains help-text "Examples:")))))

(run-tests)