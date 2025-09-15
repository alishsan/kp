(ns kp.common
  "Common mathematical functions and utilities for Kronig-Penney models.
   
   This module provides shared mathematical functions used across
   all KP model implementations (1D, 2D, extended, multilayer).")

(defn clamp
  "Clamp value x to range [min-val, max-val]."
  [x min-val max-val]
  (-> x
      (max min-val)
      (min max-val)))

(defn safe
  "Safely handle numerical edge cases by clamping to small positive values.
   Prevents division by zero and numerical overflow in trigonometric functions."
  [x]
  (let [eps 1.0e-10]
    (cond
      (Double/isNaN x) eps
      (Double/isInfinite x) (if (pos? x) 1.0e6 -1.0e6)
      (<= (Math/abs x) eps) eps
      :else x)))

(defn alpha
  "Calculate alpha = sqrt(2*m*E/ħ²) for free particle.
   In our units: alpha = sqrt(E/mu) where mu = ħ²/(2m)."
  [E mu]
  (Math/sqrt (/ E mu)))

(defn beta
  "Calculate beta = sqrt(2*m*(V0-E)/ħ²) for barrier region.
   In our units: beta = sqrt((V0-E)/mu) where mu = ħ²/(2m)."
  [V0 E mu]
  (Math/sqrt (/ (- V0 E) mu)))

(defn gamma
  "Calculate gamma = sqrt(2*m*(E-V0)/ħ²) for above-barrier region.
   In our units: gamma = sqrt((E-V0)/mu) where mu = ħ²/(2m)."
  [E V0 mu]
  (Math/sqrt (/ (- E V0) mu)))

(defn bisect
  "Bisection method for finding root of function f in interval [a, b].
   Returns approximate root within tolerance tol after max-iter iterations."
  [f a b tol max-iter]
  (let [a (double a)
        b (double b)
        tol (double tol)
        max-iter (long max-iter)]
    (loop [left a
           right b
           iter 0]
      (if (or (>= iter max-iter) (< (- right left) tol))
        (/ (+ left right) 2.0)
        (let [mid (/ (+ left right) 2.0)
              f-mid (f mid)]
          (if (<= (* (f left) f-mid) 0.0)
            (recur left mid (inc iter))
            (recur mid right (inc iter))))))))

(defn mat-mul
  "Multiply two 2x2 matrices represented as [m11 m12 m21 m22].
   Returns [m11 m12 m21 m22] of product matrix."
  [A B]
  (let [[a11 a12 a21 a22] A
        [b11 b12 b21 b22] B]
    [(+ (* a11 b11) (* a12 b21))  ; m11
     (+ (* a11 b12) (* a12 b22))  ; m12
     (+ (* a21 b11) (* a22 b21))  ; m21
     (+ (* a21 b12) (* a22 b22))])) ; m22

(defn layer-matrix
  "Calculate transfer matrix for a single layer in the multilayer model.
  
  For a layer with potential V and width w, the transfer matrix is:
  M = [cos(k*w)     sin(k*w)/k    ]
      [-k*sin(k*w)  cos(k*w)      ]
  
  where k = sqrt(2*m*(E-V)/ħ²) = sqrt((E-V)/mu)
  
  Args:
  - E: energy
  - V: potential in layer
  - w: layer width
  - mu: ħ²/(2m)
  
  Returns 2x2 transfer matrix as [m11 m12 m21 m22]"
  [E V w mu]
  (let [E (double E)
        V (double V)
        w (double w)
        mu (double mu)
        k (Math/sqrt (/ (- E V) mu))
        k (safe k)
        cos-kw (Math/cos (* k w))
        sin-kw (Math/sin (* k w))]
    [cos-kw                    ; m11
     (/ sin-kw k)             ; m12
     (* (- k) sin-kw)         ; m21
     cos-kw]))                ; m22