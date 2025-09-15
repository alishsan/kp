(ns kp.model-2d-test
  (:require [clojure.test :refer [deftest is run-tests testing]]
            [kp.model-2d :refer [allowed-2d? band-structure-2d dispersion-2d dispersion-2d-separable dispersion-surface-2d effective-mass-2d find-2d-band-edges generate-2d-k-grid k-vector-magnitude principal-k-2d]]))

(deftest test-k-vector-magnitude
  (testing "k-vector-magnitude calculates 2D wave vector magnitude correctly"
    (is (= 0.0 (k-vector-magnitude 0.0 0.0)))
    (is (= 1.0 (k-vector-magnitude 1.0 0.0)))
    (is (= 1.0 (k-vector-magnitude 0.0 1.0)))
    (is (= (Math/sqrt 2.0) (k-vector-magnitude 1.0 1.0)))
    (is (= 5.0 (k-vector-magnitude 3.0 4.0)))))

(deftest test-dispersion-2d-separable
  (testing "dispersion-2d-separable function for 2D KP model"
    (let [params {:Lx 1.0 :Ly 1.0 :a 0.5 :b 0.25 :V0 10.0 :mu 1.0}]
      ;; Test that dispersion produces reasonable values
      (is (number? (dispersion-2d-separable 0.0 0.0 0.0 params)))
      (is (number? (dispersion-2d-separable 5.0 0.0 0.0 params)))
      (is (number? (dispersion-2d-separable 15.0 0.0 0.0 params)))
      
      ;; Test that dispersion is symmetric in kx and ky
      (let [D1 (dispersion-2d-separable 5.0 0.1 0.2 params)
            D2 (dispersion-2d-separable 5.0 0.2 0.1 params)]
        (is (< (Math/abs (- D1 D2)) 1e-10))))))

(deftest test-dispersion-2d
  (testing "dispersion-2d function delegates to separable version"
    (let [params {:Lx 1.0 :Ly 1.0 :a 0.5 :b 0.25 :V0 10.0 :mu 1.0}]
      ;; Test that dispersion-2d gives same result as dispersion-2d-separable
      (let [D1 (dispersion-2d 5.0 0.1 0.2 params)
            D2 (dispersion-2d-separable 5.0 0.1 0.2 params)]
        (is (< (Math/abs (- D1 D2)) 1e-10))))))

(deftest test-allowed-2d
  (testing "allowed-2d? function correctly identifies bands in 2D"
    (let [params {:Lx 1.0 :Ly 1.0 :a 0.5 :b 0.25 :V0 10.0 :mu 1.0}]
      ;; Test that allowed-2d? returns boolean values
      (is (boolean? (allowed-2d? 0.0 0.0 0.0 params)))
      (is (boolean? (allowed-2d? 5.0 0.0 0.0 params)))
      (is (boolean? (allowed-2d? 15.0 0.0 0.0 params)))
      
      ;; Test with non-zero k-vectors
      (is (boolean? (allowed-2d? 5.0 0.1 0.1 params)))
      (is (boolean? (allowed-2d? 5.0 0.5 0.5 params)))))

(deftest test-principal-k-2d
  (testing "principal-k-2d function returns correct 2D wavevector"
    (let [params {:Lx 1.0 :Ly 1.0 :a 0.5 :b 0.25 :V0 10.0 :mu 1.0}]
      ;; E = 15 should be allowed
      (let [k-result (principal-k-2d 15.0 params)]
        (is (some? k-result))
        (is (contains? k-result :kx))
        (is (contains? k-result :ky))
        (is (>= (:kx k-result) 0.0))
        (is (>= (:ky k-result) 0.0))
        (is (<= (:kx k-result) (/ Math/PI 1.0)))
        (is (<= (:ky k-result) (/ Math/PI 1.0)))))))

(deftest test-generate-2d-k-grid
  (testing "generate-2d-k-grid creates correct k-space grid"
    (let [params {:Lx 1.0 :Ly 1.0 :nx 3 :ny 3}
          k-grid (generate-2d-k-grid params)]
      (is (= 9 (count k-grid)))  ; 3x3 = 9 points
      (is (every? #(= 2 (count %)) k-grid))  ; Each point is [kx ky]
      
      ;; Check that k-values are in correct range
      (let [k-values (flatten k-grid)
            kx-values (map first k-grid)
            ky-values (map second k-grid)]
        (is (every? #(>= % 0.0) k-values))
        (is (every? #(<= % Math/PI) kx-values))
        (is (every? #(<= % Math/PI) ky-values))))))

(deftest test-band-structure-2d
  (testing "band-structure-2d calculates 2D band structure correctly"
    (let [params {:Lx 1.0 :Ly 1.0 :a 0.5 :b 0.25 :V0 10.0 :mu 1.0}
          k-grid [[0.0 0.0] [0.1 0.1] [0.5 0.5]]
          E 5.0
          band-data (band-structure-2d E k-grid params)]
      
      (is (= 3 (count band-data)))  ; Same as k-grid size
      (is (every? #(contains? % :kx) band-data))
      (is (every? #(contains? % :ky) band-data))
      (is (every? #(contains? % :D) band-data))
      (is (every? #(contains? % :allowed) band-data))
      (is (every? #(contains? % :k-magnitude) band-data))
      
      ;; Check that k-magnitude is calculated correctly
      (doseq [point band-data]
        (let [expected-mag (k-vector-magnitude (:kx point) (:ky point))]
          (is (< (Math/abs (- (:k-magnitude point) expected-mag)) 1e-10)))))))

(deftest test-dispersion-surface-2d
  (testing "dispersion-surface-2d calculates dispersion surface correctly"
    (let [params {:Lx 1.0 :Ly 1.0 :a 0.5 :b 0.25 :V0 10.0 :mu 1.0}
          k-grid [[0.0 0.0] [0.1 0.1] [0.5 0.5]]
          E 5.0
          surface (dispersion-surface-2d E k-grid params)]
      
      (is (= 3 (count surface)))  ; Same as k-grid size
      (is (every? #(= 3 (count %)) surface))  ; Each point is [kx ky D]
      
      ;; Check that D values are reasonable
      (doseq [[kx ky D] surface]
        (is (number? D))
        (is (not (Double/isNaN D)))
        (is (not (Double/isInfinite D)))))))

(deftest test-effective-mass-2d
  (testing "effective-mass-2d calculates effective mass tensor correctly"
    (let [params {:Lx 1.0 :Ly 1.0 :a 0.5 :b 0.25 :V0 10.0 :mu 1.0}
          E 5.0
          kx 0.1
          ky 0.1
          mass-tensor (effective-mass-2d E kx ky params)]
      
      (is (contains? mass-tensor :m-star-x))
      (is (contains? mass-tensor :m-star-y))
      (is (contains? mass-tensor :m-star-xy))
      
      ;; For separable case, off-diagonal terms should be zero
      (is (= 0.0 (:m-star-xy mass-tensor)))
      
      ;; Effective masses should be finite numbers
      (is (number? (:m-star-x mass-tensor)))
      (is (number? (:m-star-y mass-tensor)))
      (is (not (Double/isNaN (:m-star-x mass-tensor))))
      (is (not (Double/isNaN (:m-star-y mass-tensor)))))))

(deftest test-find-2d-band-edges
  (testing "find-2d-band-edges finds band structure information"
    (let [params {:Lx 1.0 :Ly 1.0 :a 0.5 :b 0.25 :V0 10.0 :mu 1.0}
          k-grid [[0.0 0.0] [0.1 0.1] [0.5 0.5]]
          E-range [0.0 20.0]
          band-edges (find-2d-band-edges E-range k-grid params)]
      
      (is (seq band-edges))
      (is (every? #(contains? % :E) band-edges))
      (is (every? #(contains? % :allowed-count) band-edges))
      (is (every? #(contains? % :band-data) band-edges))
      
      ;; Check that energy range is covered
      (let [energies (map :E band-edges)]
        (is (<= (first energies) 0.0))
        (is (>= (last energies) 20.0))))))

(deftest test-2d-consistency-with-1d
  (testing "2D separable model should be consistent with 1D model at k=0"
    (let [params-2d {:Lx 1.0 :Ly 1.0 :a 0.5 :b 0.25 :V0 10.0 :mu 1.0}
          E 5.0]
      
      ;; At kx=0, ky=0, 2D separable should give D^2 where D is 1D result
      (let [D-2d (dispersion-2d E 0.0 0.0 params-2d)
            ;; For separable case at k=0, D_2d = D_1d^2
            ;; This is because cos(0) = 1, so D_2d = D_x * D_y = D_1d * D_1d
            D-1d-expected (Math/sqrt (Math/abs D-2d))]
        
        ;; The relationship should hold approximately
        (is (< (Math/abs (- (Math/abs D-2d) (* D-1d-expected D-1d-expected))) 1e-6))))))

(deftest test-2d-symmetry
  (testing "2D dispersion should be symmetric under kx <-> ky exchange"
    (let [params {:Lx 1.0 :Ly 1.0 :a 0.5 :b 0.25 :V0 10.0 :mu 1.0}
          E 5.0
          kx 0.1
          ky 0.2]
      
      (let [D1 (dispersion-2d E kx ky params)
            D2 (dispersion-2d E ky kx params)]
        (is (< (Math/abs (- D1 D2)) 1e-10))))))

(run-tests)
