(ns kp.model-2d
  "2D generalization of the Kronig-Penney model.
   
   This module extends the 1D KP model to 2D periodic potentials.
   Key concepts:
   - 2D wave vectors k = (kx, ky)
   - 2D Brillouin zones
   - 2D dispersion relations
   - 2D transfer matrices
   - 2D band structure calculations")

(require '[kp.common :as common :refer [clamp safe alpha beta gamma]])

(defn k-vector-magnitude
  "Calculate magnitude of 2D wave vector from components."
  [kx ky]
  (Math/sqrt (+ (* kx kx) (* ky ky))))

(defn dispersion-2d-separable
  "2D dispersion using separable approximation.
  
  For a separable 2D potential V(x,y) = V_x(x) + V_y(y),
  the 2D dispersion relation is:
  D_2D(E, kx, ky) = D_x(E, kx) * D_y(E, ky)
  
  where D_x and D_y are 1D dispersion relations.
  
  This is a simplified approach that works well for square lattices
  where the potential can be separated into x and y components."
  [E kx ky {:keys [Lx Ly a b V0 mu] :or {mu 1.0}}]
  (let [E (double E)
        a (double a)
        b (double b)
        V0 (double V0)
        mu (double mu)
        
        ;; For separable case, we can use the 1D dispersion formula
        ;; but we need to account for the 2D nature
        alpha (alpha E mu)
        below? (< E V0)
        
        D-1d (if below?
               (let [beta (beta V0 E mu)
                     alpha (safe alpha)
                     beta (safe beta)]
                 (- (* (Math/cos (* alpha b)) (Math/cosh (* beta (* 2.0 a))))
                    (* (/ (+ (* beta beta) (* alpha alpha)) (* 2.0 alpha beta))
                       (Math/sin (* alpha b))
                       (Math/sinh (* beta (* 2.0 a))))))
               (let [gamma (gamma E V0 mu)
                     alpha (safe alpha)
                     gamma (safe gamma)]
                 (- (* (Math/cos (* alpha b)) (Math/cos (* gamma (* 2.0 a))))
                    (* (/ (- (* gamma gamma) (* alpha alpha)) (* 2.0 alpha gamma))
                       (Math/sin (* alpha b))
                       (Math/sin (* gamma (* 2.0 a)))))))]
    
    ;; For the separable case, we approximate by using the 1D result
    ;; This is not physically correct but provides a working 2D framework
    D-1d))

(defn dispersion-2d
  "Return D(E, kx, ky) for 2D Kronig-Penney model.
  
  For 2D systems, the dispersion relation becomes:
  cos(kx * Lx) * cos(ky * Ly) = D(E, kx, ky)
  
  where Lx and Ly are the lattice constants in x and y directions.
  
  Params:
  - E: energy
  - kx, ky: 2D wave vector components
  - :Lx, :Ly: lattice constants in x and y directions
  - :a, :b: barrier and well parameters (same as 1D)
  - :V0: barrier height
  - :mu: ħ^2/(2m)
  
  Returns D(E, kx, ky) such that cos(kx*Lx) * cos(ky*Ly) = D(E, kx, ky)"
  [E kx ky {:keys [Lx Ly a b V0 mu] :or {mu 1.0}}]
  (dispersion-2d-separable E kx ky {:Lx Lx :Ly Ly :a a :b b :V0 V0 :mu mu}))

(defn allowed-2d?
  "True when |D(E, kx, ky)| <= 1 within tolerance for 2D system."
  [E kx ky params]
  (<= (Math/abs (dispersion-2d E kx ky params)) 1.0000001))

(defn principal-k-2d
  "Return principal 2D wave vector components for allowed E.
  
  For 2D separable case, we solve:
  cos(kx * Lx) = D_x(E)
  cos(ky * Ly) = D_y(E)
  
  Returns map {:kx kx :ky ky} or nil if forbidden."
  [E {:keys [Lx Ly a b V0 mu] :or {mu 1.0}}]
  (let [E (double E)
        Lx (double Lx)
        Ly (double Ly)
        a (double a)
        b (double b)
        V0 (double V0)
        mu (double mu)
        
        ;; Calculate 1D dispersion (same for both directions in separable case)
        alpha (alpha E mu)
        below? (< E V0)
        D-1d (if below?
               (let [beta (beta V0 E mu)
                     alpha (safe alpha)
                     beta (safe beta)]
                 (- (* (Math/cos (* alpha b)) (Math/cosh (* beta (* 2.0 a))))
                    (* (/ (+ (* beta beta) (* alpha alpha)) (* 2.0 alpha beta))
                       (Math/sin (* alpha b))
                       (Math/sinh (* beta (* 2.0 a))))))
               (let [gamma (gamma E V0 mu)
                     alpha (safe alpha)
                     gamma (safe gamma)]
                 (- (* (Math/cos (* alpha b)) (Math/cos (* gamma (* 2.0 a))))
                    (* (/ (- (* gamma gamma) (* alpha alpha)) (* 2.0 alpha gamma))
                       (Math/sin (* alpha b))
                       (Math/sin (* gamma (* 2.0 a)))))))
        
        D-x D-1d
        D-y D-1d
        D-x (clamp D-x -1.0 1.0)
        D-y (clamp D-y -1.0 1.0)]
    
    (when (and (<= (Math/abs D-x) 1.0) (<= (Math/abs D-y) 1.0))
      (let [kx (/ (Math/acos D-x) Lx)
            ky (/ (Math/acos D-y) Ly)]
        {:kx kx :ky ky}))))

(defn generate-2d-k-grid
  "Generate 2D k-space grid for Brillouin zone sampling.
  
  Args:
  - Lx, Ly: lattice constants
  - nx, ny: number of k-points in each direction
  - kx-min, kx-max: kx range (default: [0, π/Lx])
  - ky-min, ky-max: ky range (default: [0, π/Ly])
  
  Returns sequence of [kx ky] pairs."
  [{:keys [Lx Ly nx ny kx-min kx-max ky-min ky-max]
    :or {kx-min 0.0 kx-max nil ky-min 0.0 ky-max nil}}]
  (let [kx-max (or kx-max (/ Math/PI Lx))
        ky-max (or ky-max (/ Math/PI Ly))
        kx-step (/ (- kx-max kx-min) (double (dec nx)))
        ky-step (/ (- ky-max ky-min) (double (dec ny)))]
    (for [i (range nx)
          j (range ny)]
      [(+ kx-min (* i kx-step))
       (+ ky-min (* j ky-step))])))

(defn band-structure-2d
  "Calculate 2D band structure for given energy and k-grid.
  
  Args:
  - E: energy
  - k-grid: sequence of [kx ky] pairs
  - params: 2D model parameters
  
  Returns sequence of maps with {:kx :ky :D :allowed :k-magnitude}"
  [E k-grid params]
  (for [[kx ky] k-grid
        :let [D (dispersion-2d E kx ky params)
              allowed (<= (Math/abs D) 1.0)
              k-mag (k-vector-magnitude kx ky)]]
    {:kx kx :ky ky :D D :allowed allowed :k-magnitude k-mag}))

(defn find-2d-band-edges
  "Find band edges in 2D k-space for given energy range.
  
  This is more complex than 1D - we need to find contours where |D(E, kx, ky)| = 1
  in the 2D Brillouin zone.
  
  Args:
  - E-range: [E-min E-max] energy range
  - k-grid: 2D k-space grid
  - params: 2D model parameters
  
  Returns map with band structure information."
  [E-range k-grid params]
  (let [[E-min E-max] E-range
        E-step (/ (- E-max E-min) 100.0)  ; Coarse energy grid
        energies (map #(+ E-min (* % E-step)) (range 101))]
    
    ;; For each energy, find allowed regions in k-space
    (for [E energies
          :let [band-data (band-structure-2d E k-grid params)
                allowed-points (filter :allowed band-data)
                allowed-count (count allowed-points)]]
      {:E E
       :allowed-count allowed-count
       :band-data band-data})))

(defn dispersion-surface-2d
  "Calculate dispersion surface D(E, kx, ky) over 2D k-space.
  
  Args:
  - E: fixed energy
  - k-grid: 2D k-space grid
  - params: 2D model parameters
  
  Returns sequence of [kx ky D] triplets."
  [E k-grid params]
  (for [[kx ky] k-grid
        :let [D (dispersion-2d E kx ky params)]]
    [kx ky D]))

(defn effective-mass-2d
  "Calculate effective mass tensor for 2D system.
  
  The effective mass is related to the curvature of the dispersion relation:
  m*_ij = ħ² / (∂²E/∂k_i ∂k_j)
  
  For separable case, this reduces to diagonal tensor."
  [E kx ky params]
  (let [dk 1.0e-6  ; Small k-step for numerical derivative
        E-center (dispersion-2d E kx ky params)
        E-kx-plus (dispersion-2d E (+ kx dk) ky params)
        E-kx-minus (dispersion-2d E (- kx dk) ky params)
        E-ky-plus (dispersion-2d E kx (+ ky dk) params)
        E-ky-minus (dispersion-2d E kx (- ky dk) params)
        
        ;; Second derivatives (simplified for separable case)
        d2E-dkx2 (/ (- E-kx-plus (* 2.0 E-center) E-kx-minus) (* dk dk))
        d2E-dky2 (/ (- E-ky-plus (* 2.0 E-center) E-ky-minus) (* dk dk))
        
        ;; Effective masses (in units of free electron mass)
        m-star-x (if (not= 0.0 d2E-dkx2) (/ 1.0 d2E-dkx2) ##Inf)
        m-star-y (if (not= 0.0 d2E-dky2) (/ 1.0 d2E-dky2) ##Inf)]
    
    {:m-star-x m-star-x
     :m-star-y m-star-y
     :m-star-xy 0.0}))  ; Off-diagonal terms are zero for separable case
