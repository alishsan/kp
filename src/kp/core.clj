(ns kp.core
  (:require [clojure.tools.cli :refer [parse-opts]]
            [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [clojure.string :as str]
            [kp.model-1d :as model-1d]
            [kp.model-2d :as model-2d]))

(def cli-options
  [["-a" "--a A" "Half barrier width (barrier width = 2a)"
    :default 0.5 :parse-fn #(Double/parseDouble %)]
   ["-b" "--b B" "Well width on each side (total well width = 2b)"
    :default 0.25 :parse-fn #(Double/parseDouble %)]
   ["-V" "--V0 V0" "Barrier height V0"
    :default 10.0 :parse-fn #(Double/parseDouble %)]
   ["-m" "--mu MU" "mu = ħ^2/(2m)"
    :default 1.0 :parse-fn #(Double/parseDouble %)]
   ["-n" "--steps N" "Energy samples"
    :default 5000 :parse-fn #(Long/parseLong %)]
   ["--Emin" "--Emin EMIN" "Minimum energy"
    :default 0.0 :parse-fn #(Double/parseDouble %)]
   ["--Emax" "--Emax EMAX" "Maximum energy"
    :default 50.0 :parse-fn #(Double/parseDouble %)]
   ["--U1" "--U1 U1" "Height of first barrier (U1)"
    :default 8.0 :parse-fn #(Double/parseDouble %)]
   ["--U2" "--U2 U2" "Height of second barrier (U2, U2 >= U1)"
    :default 12.0 :parse-fn #(Double/parseDouble %)]
   ["--extended" "--extended" "Use extended KP model (U1-well-U2-well structure)"
    :default false]
   ["--layers" "--layers LAYERS" (str "Compact multilayer spec. Examples: "
                                     "'b:0.3,U:8:0.2,b:0.3' (well, barrier V=8, well) "
                                     "or 'U:12:0.2,b:0.4,U:8:0.2,b:0.4' . ")
    :default nil]
   ["--2d" "--2d" "Use 2D Kronig-Penney model"
    :default false]
   ["--Lx" "--Lx LX" "Lattice constant in x-direction (2D only)"
    :default 1.0 :parse-fn #(Double/parseDouble %)]
   ["--Ly" "--Ly LY" "Lattice constant in y-direction (2D only)"
    :default 1.0 :parse-fn #(Double/parseDouble %)]
   ["--nx" "--nx NX" "Number of k-points in x-direction (2D only)"
    :default 20 :parse-fn #(Long/parseLong %)]
   ["--ny" "--ny NY" "Number of k-points in y-direction (2D only)"
    :default 20 :parse-fn #(Long/parseLong %)]
   ["--kx-min" "--kx-min KX_MIN" "Minimum kx value (2D only)"
    :default 0.0 :parse-fn #(Double/parseDouble %)]
   ["--kx-max" "--kx-max KX_MAX" "Maximum kx value (2D only)"
    :default nil :parse-fn #(Double/parseDouble %)]
   ["--ky-min" "--ky-min KY_MIN" "Minimum ky value (2D only)"
    :default 0.0 :parse-fn #(Double/parseDouble %)]
   ["--ky-max" "--ky-max KY_MAX" "Maximum ky value (2D only)"
    :default nil :parse-fn #(Double/parseDouble %)]
   ["-o" "--out PATH" "Output CSV for dispersion samples"
    :default "dispersion.csv"]
   [nil "--help"]])

(defn usage [summary]
  (str "Kronig–Penney model CLI\n\n"
       "Options:\n" summary "\n\n"
       "Examples:\n"
       "  lein run -a 0.5 -b 0.25 -V 12.0 --Emax 40 --steps 10000 -o dispersion.csv\n"
       "  lein run -a 0.3 -b 0.2 --U1 8.0 --U2 12.0 --extended --Emax 40 -o extended.csv\n"
       "  lein run --layers b:0.4,U:12:0.2,b:0.4 --Emax 40 --steps 10000 -o multi.csv\n"
       "  lein run --2d -a 0.5 -b 0.25 -V 12.0 --Lx 1.0 --Ly 1.0 --nx 20 --ny 20 -o dispersion-2d.csv\n"))

(defn exit! [code msg]
  (binding [*out* (if (zero? code) *out* *err*)]
    (println msg))
  (System/exit code))

(defn parse-layer-token
  "Parse a single layer token: 'b:WIDTH' for V=0 or 'U:V:WIDTH' for barrier.
  Returns {:w width :V V}."
  [tok]
  (let [parts (str/split tok #":")]
    (case (first parts)
      "b" (let [[_ w] parts]
            {:w (Double/parseDouble w) :V 0.0})
      "U" (let [[_ v w] parts]
            {:w (Double/parseDouble w) :V (Double/parseDouble v)})
      ;; fallback: treat as width-only well
      {:w (Double/parseDouble (first parts)) :V 0.0})))

(defn parse-layers
  "Parse compact layers string into vector of {:w :V}."
  [s]
  (->> (str/split (str/trim s) #",+")
       (remove str/blank?)
       (map parse-layer-token)
       (vec)))

(defn generate-energy-grid
  "Generate energy grid from Emin to Emax with given number of steps."
  [Emin Emax steps]
  (let [step (/ (- Emax Emin) (double steps))]
    (map #(+ Emin (* % step)) (range (inc steps)))))

(defn write-csv-data
  "Write data rows to CSV file with given header."
  [filename header rows]
  (with-open [w (io/writer filename)]
    (csv/write-csv w (cons header rows))))

(defn process-multilayer
  "Process multilayer transfer-matrix calculation and write outputs."
  [layers energies mu out]
  (let [layers* (parse-layers layers)
        rows (for [E energies
                   :let [{:keys [D L]} (model-1d/dispersion-multilayer E layers* {:mu mu})
                         allowed (<= (Math/abs D) 1.0)
                         k (model-1d/principal-k-from-L D L)
                         kplus k
                         kminus (- k)]]
               [(format "%.12g" E)
                (format "%.12g" D)
                (if allowed "true" "false")
                (format "%.12g" kminus)
                (format "%.12g" kplus)])
        header ["E" "D" "allowed" "k_minus" "k_plus"]]
    (write-csv-data out header rows)
    (let [edges (model-1d/band-edges-multilayer layers* {:mu mu :Emin (first energies) :Emax (last energies) :steps (count energies)})
          bands-out (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")
          header2 ["E_lo" "E_hi"]
          rows2 (for [[lo hi] edges]
                  [(format "%.12g" lo) (format "%.12g" hi)])]
      (write-csv-data bands-out header2 rows2))
    (println (str "Wrote samples to " out))
    (println (str "Wrote band edges to " (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")))))

(defn process-extended-kp
  "Process extended Kronig-Penney calculation (U1-well-U2-well) and write outputs."
  [a b U1 U2 mu energies out]
  (let [params {:a a :b b :U1 U1 :U2 U2 :mu mu}
        rows (for [E energies
                   :let [{:keys [D L]} (model-1d/dispersion-extended-kp E params)
                         allowed (<= (Math/abs D) 1.0)
                         k (model-1d/principal-k-extended-kp E params)
                         kplus k
                         kminus (- k)]]
               [(format "%.12g" E)
                (format "%.12g" D)
                (if allowed "true" "false")
                (format "%.12g" kminus)
                (format "%.12g" kplus)])
        header ["E" "D" "allowed" "k_minus" "k_plus"]]
    (write-csv-data out header rows)
    (let [edges (model-1d/band-edges-extended-kp (merge params {:Emin (first energies) :Emax (last energies) :steps (count energies)}))
          bands-out (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")
          header2 ["E_lo" "E_hi"]
          rows2 (for [[lo hi] edges]
                  [(format "%.12g" lo) (format "%.12g" hi)])]
      (write-csv-data bands-out header2 rows2))
    (println (str "Wrote samples to " out))
    (println (str "Wrote band edges to " (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")))))

(defn process-2d-kp
  "Process 2D Kronig-Penney calculation and write outputs."
  [a b V0 mu Lx Ly nx ny kx-min kx-max ky-min ky-max energies out]
  (let [params {:a a :b b :V0 V0 :mu mu :Lx Lx :Ly Ly}
        kx-max (or kx-max (/ Math/PI Lx))
        ky-max (or ky-max (/ Math/PI Ly))
        k-grid-params {:Lx Lx :Ly Ly :nx nx :ny ny 
                       :kx-min kx-min :kx-max kx-max
                       :ky-min ky-min :ky-max ky-max}
        k-grid (model-2d/generate-2d-k-grid k-grid-params)
        
        ;; For 2D, we calculate band structure for each energy
        rows (for [E energies
                   [kx ky] k-grid
                   :let [D (model-2d/dispersion-2d E kx ky params)
                         allowed (<= (Math/abs D) 1.0)
                         k-mag (model-2d/k-vector-magnitude kx ky)]]
               [(format "%.12g" E)
                (format "%.12g" kx)
                (format "%.12g" ky)
                (format "%.12g" D)
                (if allowed "true" "false")
                (format "%.12g" k-mag)])
        header ["E" "kx" "ky" "D" "allowed" "k_magnitude"]]
    (write-csv-data out header rows)
    
    ;; For 2D, we also generate a dispersion surface at a fixed energy
    (let [E-fixed (/ (+ (first energies) (last energies)) 2.0)
          surface (model-2d/dispersion-surface-2d E-fixed k-grid params)
          surface-out (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-surface.csv")
          surface-header ["kx" "ky" "D"]
          surface-rows (for [[kx ky D] surface]
                        [(format "%.12g" kx)
                         (format "%.12g" ky)
                         (format "%.12g" D)])]
      (write-csv-data surface-out surface-header surface-rows))
    
    (println (str "Wrote 2D band structure to " out))
    (println (str "Wrote dispersion surface to " (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-surface.csv")))))

(defn process-single-kp
  "Process single Kronig-Penney calculation and write outputs."
  [a b V0 mu energies out]
  (let [params {:a a :b b :V0 V0 :mu mu}
        rows (for [E energies
                   :let [D (model-1d/dispersion E params)
                         allowed (<= (Math/abs D) 1.0)
                         k (model-1d/principal-k E params)
                         kplus k
                         kminus (- k)]]
               [(format "%.12g" E)
                (format "%.12g" D)
                (if allowed "true" "false")
                (format "%.12g" kminus)
                (format "%.12g" kplus)])
        header ["E" "D" "allowed" "k_minus" "k_plus"]]
    (write-csv-data out header rows)
    (let [edges (model-1d/band-edges (merge params {:Emin (first energies) :Emax (last energies) :steps (count energies)}))
          bands-out (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")
          header2 ["E_lo" "E_hi"]
          rows2 (for [[lo hi] edges]
                  [(format "%.12g" lo) (format "%.12g" hi)])]
      (write-csv-data bands-out header2 rows2))
    (println (str "Wrote samples to " out))
    (println (str "Wrote band edges to " (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")))))

(defn validate-extended-kp-params
  "Validate extended KP parameters."
  [a b U1 U2]
  (when (or (<= a 0.0) (<= b 0.0))
    (exit! 1 (str "Require a > 0 and b > 0, got a=" a ", b=" b)))
  (when (< U2 U1)
    (exit! 1 (str "Require U2 >= U1, got U1=" U1 ", U2=" U2))))

(defn validate-single-kp-params
  "Validate single KP parameters."
  [a b]
  (when (or (<= a 0.0) (<= b 0.0))
    (exit! 1 (str "Require a > 0 and b > 0, got a=" a ", b=" b))))

(defn -main [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args cli-options)
        {:keys [a b V0 mu steps Emin Emax out layers U1 U2 extended 
                Lx Ly nx ny kx-min kx-max ky-min ky-max]} options]
    (when (:help options) (exit! 0 (usage summary)))
    (when (seq errors) (exit! 1 (str (str/join "\n" errors) "\n\n" (usage summary))))
    
    (let [energies (generate-energy-grid Emin Emax steps)]
      (cond
        (:2d options)
        ;; 2D KP mode
        (do
          (validate-single-kp-params a b)
          (process-2d-kp a b V0 mu Lx Ly nx ny kx-min kx-max ky-min ky-max energies out))
        
        (some? layers)
        ;; Multilayer transfer-matrix mode
        (process-multilayer layers energies mu out)
        
        extended
        ;; Extended KP mode (U1-well-U2-well)
        (do
          (validate-extended-kp-params a b U1 U2)
          (process-extended-kp a b U1 U2 mu energies out))
        
        :else
        ;; Single KP closed-form mode
        (do
          (validate-single-kp-params a b)
          (process-single-kp a b V0 mu energies out))))))
