(ns kp.model)

(def ^:private eps 1.0e-12)

(defn clamp
  "Clamp x into [lo, hi]."
  [x lo hi]
  (-> x (max lo) (min hi)))

(defn- safe
  "Avoid exact zeros to reduce division/acos domain errors."
  [x]
  (let [ax (Math/abs x)]
    (if (< ax eps)
      (if (neg? x) (- eps) eps)
      x)))

(defn- alpha
  "Wave number in the well region (V=0); mu = ħ^2/(2m)."
  [E mu]
  (let [E (max 0.0 E)]
    (Math/sqrt (/ E (max eps mu)))))

(defn- beta
  "E < V0: evanescent wave number in barrier; mu = ħ^2/(2m)."
  [V0 E mu]
  (let [d (- V0 E)]
    (if (pos? d)
      (Math/sqrt (/ d (max eps mu)))
      0.0)))

(defn- gamma
  "E > V0: oscillatory wave number in barrier; mu = ħ^2/(2m)."
  [E V0 mu]
  (let [d (- E V0)]
    (if (pos? d)
      (Math/sqrt (/ d (max eps mu)))
      0.0)))

(defn dispersion
  "Return D(E) such that cos(k a) = D(E) for the Kronig–Penney model.

  Params map keys:
  - :a  (lattice period)
  - :b  (well width, 0 < b < a)
  - :V0 (barrier height)
  - :mu (ħ^2/(2m)); default 1.0

  For E < V0 uses cosh/sinh form; for E > V0 uses trig-trig form.
  The function is numerically stabilized with small eps clamps."
  [E {:keys [a b V0 mu] :or {mu 1.0}}]
  (let [a (double a)
        b (double b)
        V0 (double V0)
        mu (double mu)
        c (- a b)
        E (double E)
        alpha (alpha E mu)
        below? (< E V0)
        D (if below?
            (let [beta (beta V0 E mu)
                  alpha (safe alpha)
                  beta (safe beta)]
              (- (* (Math/cos (* alpha b)) (Math/cosh (* beta c)))
                 (* (/ (+ (* beta beta) (* alpha alpha)) (* 2.0 alpha beta))
                    (Math/sin (* alpha b))
                    (Math/sinh (* beta c)))))
            (let [gamma (gamma E V0 mu)
                  alpha (safe alpha)
                  gamma (safe gamma)]
              (- (* (Math/cos (* alpha b)) (Math/cos (* gamma c)))
                 (* (/ (- (* gamma gamma) (* alpha alpha)) (* 2.0 alpha gamma))
                    (Math/sin (* alpha b))
                    (Math/sin (* gamma c))))))]
    D))

(defn allowed?
  "True when |D(E)| <= 1 within tolerance."
  [E params]
  (<= (Math/abs (dispersion E params)) 1.0000001))

(defn principal-k
  "Return principal |k| in [0, π/a] for allowed E using k = acos(D)/a.
  Returns nil if forbidden (|D|>1)."
  [E {:keys [a] :as params}]
  (let [D (dispersion E params)
        D (clamp D -1.0 1.0)
        k (/ (Math/acos D) (double a))]
    k))

(defn- bisect
  "Bisection on function f between xlo and xhi where f(xlo) and f(xhi) have opposite signs.
  Returns root with absolute tolerance tol or after max-it iterations."
  [f xlo xhi tol max-it]
  (loop [lo (double xlo)
         hi (double xhi)
         it 0]
    (let [mid (+ lo (* 0.5 (- hi lo)))
          fmid (f mid)
          flo (f lo)]
      (cond
        (or (<= (Math/abs (- hi lo)) tol) (>= it max-it)) mid
        (<= (* fmid flo) 0.0) (recur lo mid (inc it))
        :else (recur mid hi (inc it))))))

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
        ;; end of grid
        (if in-band?
          (conj acc [start-E (double Emax)])
          acc)))))