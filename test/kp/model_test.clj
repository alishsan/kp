(ns kp.model-test
  (:require [clojure.test :refer :all]
            [kp.model :refer :all]))

(deftest test-clamp
  (testing "clamp function bounds values correctly"
    (is (= 5.0 (clamp 5.0 0.0 10.0)))
    (is (= 0.0 (clamp -5.0 0.0 10.0)))
    (is (= 10.0 (clamp 15.0 0.0 10.0)))))

(deftest test-dispersion-single-kp
  (testing "dispersion function for single KP model"
    (let [params {:a 1.0 :b 0.5 :V0 10.0 :mu 1.0}]
      ;; Test at E = 0 (should be in band for these parameters)
      (is (<= (Math/abs (dispersion 0.0 params)) 1.0))
      ;; Test at E = 5 (should be in band for these parameters)
      (is (<= (Math/abs (dispersion 5.0 params)) 1.0))
      ;; Test at E = 15 (above barrier, should be in band)
      (is (<= (Math/abs (dispersion 15.0 params)) 1.0)))))

(deftest test-allowed
  (testing "allowed? function correctly identifies bands"
    (let [params {:a 1.0 :b 0.5 :V0 10.0 :mu 1.0}]
      (is (allowed? 0.0 params))
      (is (allowed? 5.0 params))
      (is (allowed? 15.0 params)))))

(deftest test-principal-k
  (testing "principal-k function returns correct wavevector"
    (let [params {:a 1.0 :b 0.5 :V0 10.0 :mu 1.0}]
      ;; E = 15 should be allowed
      (let [k (principal-k 15.0 params)]
        (is (some? k))
        (is (>= k 0.0))
        (is (<= k (/ Math/PI 1.0)))))))

(deftest test-dispersion-multilayer
  (testing "dispersion-multilayer function computes D(E) correctly"
    (let [layers [{:w 0.4 :V 0.0} {:w 0.2 :V 10.0} {:w 0.4 :V 0.0}]
          opts {:mu 1.0}]
      ;; Test at E = 5 (should be in gap)
      (let [{:keys [D L]} (dispersion-multilayer 5.0 layers opts)]
        (is (number? D))
        (is (number? L))
        (is (= 1.0 L)))))) ;; 0.4 + 0.2 + 0.4

(deftest test-principal-k-from-L
  (testing "principal-k-from-L function computes k correctly"
    (let [D 0.5 L 2.0]
      (let [k (principal-k-from-L D L)]
        (is (some? k))
        (is (>= k 0.0))
        (is (<= k (/ Math/PI L)))))))

(deftest test-negative-potentials
  (testing "system handles negative potentials correctly"
    (let [layers [{:w 0.4 :V -5.0} {:w 0.2 :V 8.0} {:w 0.4 :V 0.0}]
          opts {:mu 1.0}]
      ;; Test at E = -3 (above negative well)
      (let [{:keys [D L]} (dispersion-multilayer -3.0 layers opts)]
        (is (number? D))
        (is (number? L))
        (is (= 1.0 L))))))

(deftest test-band-edges
  (testing "band-edges function finds band boundaries"
    (let [params {:a 1.0 :b 0.5 :V0 10.0 :mu 1.0}
          edges (band-edges {:a 1.0 :b 0.5 :V0 10.0 :mu 1.0 :Emin 0.0 :Emax 20.0 :steps 1000})]
      (is (vector? edges))
      (is (every? vector? edges))
      (is (every? #(= 2 (count %)) edges)))))

(deftest test-band-edges-multilayer
  (testing "band-edges-multilayer function finds band boundaries"
    (let [layers [{:w 0.4 :V 0.0} {:w 0.2 :V 10.0} {:w 0.4 :V 0.0}]
          edges (band-edges-multilayer layers {:mu 1.0 :Emin 0.0 :Emax 20.0 :steps 1000})]
      (is (vector? edges))
      (is (every? vector? edges))
      (is (every? #(= 2 (count %)) edges)))))

(run-tests)
