(ns kp.core
  (:require [clojure.tools.cli :refer [parse-opts]]
            [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [clojure.string :as str]
            [kp.model :as model]))

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
   ["--2d" "--2d" "Use 2D Kronig-Penney model"
    :default false]
   ["--ax" "--ax AX" "Period in x-direction for 2D model"
    :default 1.0 :parse-fn #(Double/parseDouble %)]
   ["--ay" "--ay AY" "Period in y-direction for 2D model"
    :default 1.0 :parse-fn #(Double/parseDouble %)]
   ["--bx" "--bx BX" "Barrier width in x-direction for 2D model (0 < bx < ax)"
    :default 0.5 :parse-fn #(Double/parseDouble %)]
   ["--by" "--by BY" "Barrier width in y-direction for 2D model (0 < by < ay)"
    :default 0.5 :parse-fn #(Double/parseDouble %)]
   ["--kx-steps" "--kx-steps N" "Number of kx points for 2D band structure"
    :default 50 :parse-fn #(Long/parseLong %)]
   ["--ky-steps" "--ky-steps N" "Number of ky points for 2D band structure"
    :default 50 :parse-fn #(Long/parseLong %)]
   ["--layers" "--layers LAYERS" (str "Compact multilayer spec. Examples: "
                                     "'b:0.3,U:8:0.2,b:0.3' (well, barrier V=8, well) "
                                     "or 'U:12:0.2,b:0.4,U:8:0.2,b:0.4' . ")
    :default nil]
   ["-o" "--out PATH" "Output CSV for dispersion samples"
    :default "dispersion.csv"]
   [nil "--help"]])

(defn usage [summary]
  (str "Kronig–Penney model CLI\n\n"
       "Options:\n" summary "\n\n"
       "Examples:\n"
       "  lein run -a 0.5 -b 0.25 -V 12.0 --Emax 40 --steps 10000 -o dispersion.csv\n"
       "  lein run -a 0.3 -b 0.2 --U1 8.0 --U2 12.0 --extended --Emax 40 -o extended.csv\n"
       "  lein run --ax 1.0 --ay 1.0 --bx 0.5 --by 0.5 -V 10.0 --2d --kx-steps 50 --ky-steps 50 -o 2d.csv\n"
       "  lein run --layers b:0.4,U:12:0.2,b:0.4 --Emax 40 --steps 10000 -o multi.csv\n"))

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
                   :let [{:keys [D L]} (model/dispersion-multilayer E layers* {:mu mu})
                         allowed (<= (Math/abs D) 1.0)
                         k (model/principal-k-from-L D L)
                         kplus k
                         kminus (- k)]]
               [(format "%.12g" E)
                (format "%.12g" D)
                (if allowed "true" "false")
                (format "%.12g" kminus)
                (format "%.12g" kplus)])
        header ["E" "D" "allowed" "k_minus" "k_plus"]]
    (write-csv-data out header rows)
    (let [edges (model/band-edges-multilayer layers* {:mu mu :Emin (first energies) :Emax (last energies) :steps (count energies)})
          bands-out (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")
          header2 ["E_lo" "E_hi"]
          rows2 (for [[lo hi] edges]
                  [(format "%.12g" lo) (format "%.12g" hi)])]
      (write-csv-data bands-out header2 rows2))
    (println (str "Wrote samples to " out))
    (println (str "Wrote band edges to " (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")))))

(defn process-2d-kp
  "Process 2D Kronig-Penney calculation and write outputs."
  [ax ay bx by V0 mu kx-steps ky-steps out]
  (let [params {:ax ax :ay ay :bx bx :by by :V0 V0 :mu mu}
        ;; Calculate k-space ranges (first Brillouin zone)
        kx-min (- (/ Math/PI ax))
        kx-max (/ Math/PI ax)
        ky-min (- (/ Math/PI ay))
        ky-max (/ Math/PI ay)
        kx-range [kx-min kx-max kx-steps]
        ky-range [ky-min ky-max ky-steps]
        ;; Calculate 2D band structure
        band-data (model/band-structure-2d kx-range ky-range params)
        rows (for [[kx ky E in-bz] band-data]
               [(format "%.12g" kx)
                (format "%.12g" ky)
                (format "%.12g" E)
                (if in-bz "true" "false")])
        header ["kx" "ky" "E" "in_brillouin_zone"]]
    (write-csv-data out header rows)
    (println (str "Wrote 2D band structure to " out))
    (println (str "Grid: " kx-steps " x " ky-steps " points"))
    (println (str "k-space: kx ∈ [" kx-min ", " kx-max "], ky ∈ [" ky-min ", " ky-max "]"))))

(defn process-extended-kp
  "Process extended Kronig-Penney calculation (U1-well-U2-well) and write outputs."
  [a b U1 U2 mu energies out]
  (let [params {:a a :b b :U1 U1 :U2 U2 :mu mu}
        rows (for [E energies
                   :let [{:keys [D L]} (model/dispersion-extended-kp E params)
                         allowed (<= (Math/abs D) 1.0)
                         k (model/principal-k-extended-kp E params)
                         kplus k
                         kminus (- k)]]
               [(format "%.12g" E)
                (format "%.12g" D)
                (if allowed "true" "false")
                (format "%.12g" kminus)
                (format "%.12g" kplus)])
        header ["E" "D" "allowed" "k_minus" "k_plus"]]
    (write-csv-data out header rows)
    (let [edges (model/band-edges-extended-kp (merge params {:Emin (first energies) :Emax (last energies) :steps (count energies)}))
          bands-out (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")
          header2 ["E_lo" "E_hi"]
          rows2 (for [[lo hi] edges]
                  [(format "%.12g" lo) (format "%.12g" hi)])]
      (write-csv-data bands-out header2 rows2))
    (println (str "Wrote samples to " out))
    (println (str "Wrote band edges to " (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")))))

(defn process-single-kp
  "Process single Kronig-Penney calculation and write outputs."
  [a b V0 mu energies out]
  (let [params {:a a :b b :V0 V0 :mu mu}
        rows (for [E energies
                   :let [D (model/dispersion E params)
                         allowed (<= (Math/abs D) 1.0)
                         k (model/principal-k E params)
                         kplus k
                         kminus (- k)]]
               [(format "%.12g" E)
                (format "%.12g" D)
                (if allowed "true" "false")
                (format "%.12g" kminus)
                (format "%.12g" kplus)])
        header ["E" "D" "allowed" "k_minus" "k_plus"]]
    (write-csv-data out header rows)
    (let [edges (model/band-edges (merge params {:Emin (first energies) :Emax (last energies) :steps (count energies)}))
          bands-out (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")
          header2 ["E_lo" "E_hi"]
          rows2 (for [[lo hi] edges]
                  [(format "%.12g" lo) (format "%.12g" hi)])]
      (write-csv-data bands-out header2 rows2))
    (println (str "Wrote samples to " out))
    (println (str "Wrote band edges to " (str (if (.endsWith out ".csv") (subs out 0 (- (count out) 4)) out) "-bands.csv")))))

(defn validate-2d-kp-params
  "Validate 2D KP parameters."
  [ax ay bx by]
  (when (or (<= ax 0.0) (<= ay 0.0))
    (exit! 1 (str "Require ax > 0 and ay > 0, got ax=" ax ", ay=" ay)))
  (when (or (<= bx 0.0) (<= by 0.0))
    (exit! 1 (str "Require bx > 0 and by > 0, got bx=" bx ", by=" by)))
  (when (>= bx ax)
    (exit! 1 (str "Require bx < ax, got ax=" ax ", bx=" bx)))
  (when (>= by ay)
    (exit! 1 (str "Require by < ay, got ay=" ay ", by=" by))))

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
        {:keys [a b V0 mu steps Emin Emax out layers U1 U2 extended ax ay bx by kx-steps ky-steps]} options]
    (when (:help options) (exit! 0 (usage summary)))
    (when (seq errors) (exit! 1 (str (str/join "\n" errors) "\n\n" (usage summary))))
    
    (cond
      (:2d options)
      ;; 2D KP mode
      (do
        (validate-2d-kp-params ax ay bx by)
        (process-2d-kp ax ay bx by V0 mu kx-steps ky-steps out))
      
      (some? layers)
      ;; Multilayer transfer-matrix mode
      (let [energies (generate-energy-grid Emin Emax steps)]
        (process-multilayer layers energies mu out))
      
      extended
      ;; Extended KP mode (U1-well-U2-well)
      (let [energies (generate-energy-grid Emin Emax steps)]
        (validate-extended-kp-params a b U1 U2)
        (process-extended-kp a b U1 U2 mu energies out))
      
      :else
      ;; Single KP closed-form mode
      (let [energies (generate-energy-grid Emin Emax steps)]
        (validate-single-kp-params a b)
        (process-single-kp a b V0 mu energies out)))))
