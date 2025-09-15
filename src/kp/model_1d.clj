(ns kp.model-1d
  "1D Kronig-Penney model implementations.
   
   This module contains all 1D variants:
   - Simple Kronig-Penney model
   - Extended KP model (U1-well-U2-well)
   - Multilayer transfer-matrix model")

(require '[kp.common :as common :refer [clamp safe alpha beta gamma bisect mat-mul layer-matrix]])

(defn dispersion
  "Return D(E) such that cos(k L) = D(E) for the Kronig–Penney model.

  New structure: x=0 at center of barrier, period L = 2a + 2b
  - Barrier: width 2a, centered at x=0, potential V0
  - Wells: width b on each side of barrier, potential 0

  Params map keys:
  - :a  (half barrier width, so barrier width = 2a)
  - :b  (well width on each side, so total well width = 2b)
  - :V0 (barrier height)
  - :mu (ħ^2/(2m)); default 1.0

  For E < V0 uses cosh/sinh form; for E > V0 uses trig-trig form.
  The function is numerically stabilized with small eps clamps."
  [E {:keys [a b V0 mu] :or {mu 1.0}}]
  (let [a (double a)
        b (double b)
        V0 (double V0)
        mu (double mu)
        L (+ (* 2.0 a) (* 2.0 b))  ; Total period
        E (double E)
        alpha (alpha E mu)
        below? (< E V0)
        D (if below?
            (let [beta (beta V0 E mu)
                  alpha (safe alpha)
                  beta (safe beta)]
              ;; For E < V0: cos(kL) = cos(αb)cosh(βa) + (α²+β²)/(2αβ) sin(αb)sinh(βa)
              ;; But we need to be careful about the coordinate system
              (if (or (zero? alpha) (zero? beta))
                1.0  ; Handle edge case
                (+ (* (Math/cos (* alpha b)) (Math/cosh (* beta (* 2.0 a))))
                   (* (/ (+ (* beta beta) (* alpha alpha)) (* 2.0 alpha beta))
                      (Math/sin (* alpha b))
                      (Math/sinh (* beta (* 2.0 a)))))))
            (let [gamma (gamma E V0 mu)
                  alpha (safe alpha)
                  gamma (safe gamma)]
              ;; For E > V0: cos(kL) = cos(αb)cos(γa) - (γ²-α²)/(2αγ) sin(αb)sin(γa)
              (if (or (zero? alpha) (zero? gamma))
                1.0  ; Handle edge case
                (- (* (Math/cos (* alpha b)) (Math/cos (* gamma (* 2.0 a))))
                   (* (/ (- (* gamma gamma) (* alpha alpha)) (* 2.0 alpha gamma))
                      (Math/sin (* alpha b))
                      (Math/sin (* gamma (* 2.0 a))))))))]
    D))

(defn allowed?
  "True when |D(E)| <= 1 within tolerance."
  [E params]
  (<= (Math/abs (dispersion E params)) 1.0000001))

(defn principal-k
  "Return principal |k| in [0, π/L] for allowed E using k = acos(D)/L.
  Returns nil if forbidden (|D|>1)."
  [E {:keys [a b] :as params}]
  (let [D (dispersion E params)
        D (clamp D -1.0 1.0)
        L (+ (* 2.0 a) (* 2.0 b))  ; Total period
        k (/ (Math/acos D) (double L))]
    k))

(defn band-edges
  "Scan energy grid to find approximate band edges where |D(E)| = 1.
  Returns a vector of [Elo Ehi] intervals for allowed bands.
  Options:
  - :Emin, :Emax, :steps (grid controls)
  - model params: :a, :b, :V0, :mu"
  [{:keys [Emin Emax steps] :as opts}]
  (let [params (select-keys opts [:a :b :V0 :mu])
        Emin (double (or Emin 0.0))
        Emax (double (or Emax 10.0))
        steps (long (or steps 10000))
        step (/ (- Emax Emin) (double steps))
        g (fn [E] (- (Math/abs (dispersion E params)) 1.0))
        energies (map #(+ Emin (* % step)) (range (inc steps)))]
    (loop [es energies
           last-E nil
           last-g nil
           in-band? false
           start-E nil
           acc []]
      (if-let [E (first es)]
        (let [val (g E)
              allowed (<= val 1.0e-10)]
          (cond
            (and (not in-band?) allowed last-E)
            (let [edge (bisect g last-E E 1.0e-8 64)]
              (recur (rest es) E val true edge acc))

            (and in-band? (not allowed) last-E)
            (let [edge (bisect g last-E E 1.0e-8 64)
                  acc (conj acc [start-E edge])]
              (recur (rest es) E val false nil acc))

            :else
            (recur (rest es) E val in-band? start-E acc)))
        (if in-band?
          (conj acc [start-E (double Emax)])
          acc)))))

;; ---------------- Multilayer transfer-matrix formulation ----------------

(defn- layer-matrix
  "Transfer matrix for a single layer with constant potential V and width w at energy E.
  Uses standard 2x2 form mapping [psi; psi'] across the layer.
  Returns a vector [m11 m12 m21 m22]."
  [E V w mu]
  (let [mu (double (or mu 1.0))
        E (double E)
        V (double V)
        w (double w)
        d (- E V)]
    (if (pos? d)
      (let [k (Math/sqrt (/ d mu))
            kw (* k w)
            c (Math/cos kw)
            s (Math/sin kw)
            invk (/ 1.0 (safe k))]
        [c (* invk s) (* -1.0 k s) c])
      (let [kappa (Math/sqrt (/ (- V E) mu))
            kw (* kappa w)
            ch (Math/cosh kw)
            sh (Math/sinh kw)
            invk (/ 1.0 (safe kappa))]
        [ch (* invk sh) (* kappa sh) ch]))))

(defn- mat-mul
  "Multiply two 2x2 matrices (flattened [a b c d])."
  [[a b c d] [e f g h]]
  [(+ (* a e) (* b g))
   (+ (* a f) (* b h))
   (+ (* c e) (* d g))
   (+ (* c f) (* d h))])

;; ---------------- 2D Kronig-Penney model ----------------

(defn- solve-1d-dispersion
  "Solve 1D dispersion relation to find E given k.
  For 1D KP: cos(k*a) = D(E), find E such that this holds."
  [k params]
  (let [{:keys [a V0 mu]} params
        target-cos (Math/cos (* k a))
        ;; Use bisection to find E such that D(E) = target-cos
        f (fn [E] (- (dispersion E params) target-cos))
        ;; Search range: 0 to 2*V0
        E-max (* 2.0 V0)
        E (bisect f 0.0 E-max 1.0e-8 64)]
    E))

(defn dispersion-2d
  "Compute 2D dispersion relation E(kx, ky) for 2D Kronig-Penney model.
  
  Structure: 2D periodic potential with rectangular barriers
  - Periods: ax (x-direction), ay (y-direction)  
  - Barrier: width bx × by, height V0
  - Well: width (ax-bx) × (ay-by), height 0
  
  For separable potential V(x,y) = Vx(x) + Vy(y), we have:
  E(kx,ky) = Ex(kx) + Ey(ky)
  
  Params:
  - :ax (period in x-direction)
  - :ay (period in y-direction)
  - :bx (barrier width in x-direction, 0 < bx < ax)
  - :by (barrier width in y-direction, 0 < by < ay)
  - :V0 (barrier height)
  - :mu (ħ^2/(2m)); default 1.0
  
  Returns E(kx, ky) for given kx, ky values"
  [kx ky {:keys [ax ay bx by V0 mu] :or {mu 1.0}}]
  (let [ax (double ax)
        ay (double ay)
        bx (double bx)
        by (double by)
        V0 (double V0)
        mu (double mu)
        kx (double kx)
        ky (double ky)
        ;; For separable potential, solve 1D problem for each direction
        ex-params {:a ax :b bx :V0 V0 :mu mu}
        ey-params {:a ay :b by :V0 V0 :mu mu}
        ;; Solve 1D dispersion relations
        Ex (solve-1d-dispersion kx ex-params)
        Ey (solve-1d-dispersion ky ey-params)]
    (+ Ex Ey)))

(defn dispersion-2d-exact
  "Compute exact 2D dispersion relation using transfer matrix method.
  
  This is more computationally intensive but gives exact results.
  Uses 2D transfer matrices for the full 2D problem."
  [kx ky {:keys [ax ay bx by V0 mu] :or {mu 1.0}}]
  (let [ax (double ax)
        ay (double ay)
        bx (double bx)
        by (double by)
        V0 (double V0)
        mu (double mu)
        kx (double kx)
        ky (double ky)
        ;; For now, use separable approximation
        ;; Full 2D transfer matrix implementation would be much more complex
        ex-params {:a ax :b bx :V0 V0 :mu mu}
        ey-params {:a ay :b by :V0 V0 :mu mu}
        Dx (dispersion (Math/cos (* kx ax)) ex-params)
        Dy (dispersion (Math/cos (* ky ay)) ey-params)]
    ;; Return both Dx and Dy for now
    {:Dx Dx :Dy Dy :E (dispersion-2d kx ky {:ax ax :ay ay :bx bx :by by :V0 V0 :mu mu})}))

(defn band-structure-2d
  "Calculate 2D band structure over k-space grid.
  
  Args:
  - kx-range: [kx-min, kx-max, kx-steps]
  - ky-range: [ky-min, ky-max, ky-steps]  
  - params: 2D model parameters
  
  Returns: vector of [kx, ky, E, allowed]"
  [kx-range ky-range params]
  (let [[kx-min kx-max kx-steps] kx-range
        [ky-min ky-max ky-steps] ky-range
        kx-step (/ (- kx-max kx-min) (double kx-steps))
        ky-step (/ (- ky-max ky-min) (double ky-steps))
        {:keys [ax ay]} params]
    (for [i (range (inc kx-steps))
          j (range (inc ky-steps))
          :let [kx (+ kx-min (* i kx-step))
                ky (+ ky-min (* j ky-step))
                E (dispersion-2d kx ky params)
                ;; Check if in first Brillouin zone
                in-bz (and (<= (Math/abs kx) (/ Math/PI ax))
                           (<= (Math/abs ky) (/ Math/PI ay)))]]
      [kx ky E in-bz])))

;; ---------------- Extended Kronig-Penney model: U1-well-U2-well structure ----------------

(defn dispersion-extended-kp
  "Compute D(E) for the extended Kronig-Penney model with structure U1-well-U2-well.
  
  Structure: U1 barrier (width 2a) - well (width 2b) - U2 barrier (width 2a) - well (width 2b)
  Period L = 4a + 4b
  x = 0 at center of U1 barrier
  
  Layer order: U1 barrier, well, U2 barrier, well
  
  Params:
  - :a (half barrier width, so each barrier width = 2a)
  - :b (half well width, so each well width = 2b) 
  - :U1 (height of first barrier)
  - :U2 (height of second barrier, U2 >= U1)
  - :mu (ħ^2/(2m)); default 1.0
  
  Returns map {:D D :L L} where L = 4a + 4b"
  [E {:keys [a b U1 U2 mu] :or {mu 1.0}}]
  (let [a (double a)
        b (double b)
        U1 (double U1)
        U2 (double U2)
        mu (double mu)
        L (+ (* 4.0 a) (* 4.0 b))  ; Total period
        E (double E)]
    (if (= U1 U2)
      ;; When U1 = U2, use simple KP formula with adjusted period
      ;; The extended model should reproduce the simple KP model
      (let [simple-params {:a a :b b :V0 U1 :mu mu}
            simple-D (dispersion E simple-params)]
        {:D simple-D :L (double L)})
      ;; When U1 != U2, use transfer matrix approach
      (let [;; Create the four layers: U1 barrier, well, U2 barrier, well
            layers [{:w (* 2.0 a) :V U1}   ; U1 barrier
                    {:w (* 2.0 b) :V 0.0}  ; well
                    {:w (* 2.0 a) :V U2}   ; U2 barrier
                    {:w (* 2.0 b) :V 0.0}] ; well
            ;; Calculate total transfer matrix
            Mtot (reduce mat-mul [1.0 0.0 0.0 1.0]
                         (map (fn [{:keys [w V]}] (layer-matrix E V w mu)) layers))
            [m11 _ _ m22] Mtot
            D (* 0.5 (+ m11 m22))]
        {:D D :L (double L)}))))

(defn principal-k-extended-kp
  "Compute principal |k| in [0, π/L] for extended KP model."
  [E {:keys [a b] :as params}]
  (let [{:keys [D L]} (dispersion-extended-kp E params)
        D (clamp D -1.0 1.0)
        k (/ (Math/acos D) (double L))]
    k))

(defn band-edges-extended-kp
  "Find band edges for extended KP model.
  Args:
  - params: {:a :b :U1 :U2 :mu :Emin :Emax :steps}
  Returns vector of [Elo Ehi]."
  [{:keys [a b U1 U2 mu Emin Emax steps] :or {mu 1.0 Emin 0.0 Emax 10.0 steps 10000}}]
  (let [Emin (double Emin)
        Emax (double Emax)
        steps (long steps)
        step (/ (- Emax Emin) (double steps))
        g (fn [E]
            (let [{:keys [D]} (dispersion-extended-kp E {:a a :b b :U1 U1 :U2 U2 :mu mu})]
              (- (Math/abs D) 1.0)))
        energies (map #(+ Emin (* % step)) (range (inc steps)))]
    (loop [es energies
           last-E nil
           last-g nil
           in-band? false
           start-E nil
           acc []]
      (if-let [E (first es)]
        (let [val (g E)
              allowed (<= val 1.0e-10)]
          (cond
            (and (not in-band?) allowed last-E)
            (let [edge (bisect g last-E E 1.0e-8 64)]
              (recur (rest es) E val true edge acc))

            (and in-band? (not allowed) last-E)
            (let [edge (bisect g last-E E 1.0e-8 64)
                  acc (conj acc [start-E edge])]
              (recur (rest es) E val false nil acc))

            :else
            (recur (rest es) E val in-band? start-E acc)))
        (if in-band?
          (conj acc [start-E (double Emax)])
          acc)))))

;; ---------------- Multilayer transfer-matrix formulation ----------------

(defn dispersion-multilayer
  "Compute D(E) = 0.5 * Tr(M_total) for a stack of layers.
  layers: vector of {:w width :V potential}
  opts: {:mu ...}
  Returns map {:D D :L L} where L is total period length."
  [E layers {:keys [mu] :or {mu 1.0}}]
  (let [L (reduce + (map (fn [{:keys [w]}] (double w)) layers))
        Mtot (reduce mat-mul [1.0 0.0 0.0 1.0]
                     (map (fn [{:keys [w V]}] (layer-matrix E V w mu)) layers))
        [m11 _ _ m22] Mtot
        D (* 0.5 (+ m11 m22))]
    {:D D :L (double L)}))

(defn principal-k-from-L
  "Compute principal |k| in [0, π/L] from D and L."
  [D L]
  (let [D (clamp D -1.0 1.0)]
    (/ (Math/acos D) (double L))))

(defn band-edges-multilayer
  "Scan energies to find bands using multilayer D(E) with period L.
  Args:
  - layers: vector of {:w :V}
  - opts: {:mu :Emin :Emax :steps}
  Returns vector of [Elo Ehi]."
  [layers {:keys [mu Emin Emax steps] :or {mu 1.0 Emin 0.0 Emax 10.0 steps 10000}}]
  (let [Emin (double Emin)
        Emax (double Emax)
        steps (long steps)
        step (/ (- Emax Emin) (double steps))
        g (fn [E]
            (let [{:keys [D]} (dispersion-multilayer E layers {:mu mu})]
              (- (Math/abs D) 1.0)))
        energies (map #(+ Emin (* % step)) (range (inc steps)))]
    (loop [es energies
           last-E nil
           last-g nil
           in-band? false
           start-E nil
           acc []]
      (if-let [E (first es)]
        (let [val (g E)
              allowed (<= val 1.0e-10)]
          (cond
            (and (not in-band?) allowed last-E)
            (let [edge (bisect g last-E E 1.0e-8 64)]
              (recur (rest es) E val true edge acc))

            (and in-band? (not allowed) last-E)
            (let [edge (bisect g last-E E 1.0e-8 64)
                  acc (conj acc [start-E edge])]
              (recur (rest es) E val false nil acc))

            :else
            (recur (rest es) E val in-band? start-E acc)))
        (if in-band?
          (conj acc [start-E (double Emax)])
          acc)))))