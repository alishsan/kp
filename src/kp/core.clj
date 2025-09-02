(ns kp.core
  (:require [clojure.tools.cli :refer [parse-opts]]
            [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [clojure.string :as str]
            [kp.model :as model]))

(def cli-options
  [["-a" "--a A" "Lattice period a"
    :default 1.0 :parse-fn #(Double/parseDouble %)]
   ["-b" "--b B" "Well width b (0 < b < a)"
    :default 0.5 :parse-fn #(Double/parseDouble %)]
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
       "  lein run -a 1.0 -b 0.4 -V 12.0 --Emax 40 --steps 10000 -o dispersion.csv\n"
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

(defn validate-single-kp-params
  "Validate single KP parameters."
  [a b]
  (when (not (< 0.0 b a))
    (exit! 1 (str "Require 0 < b < a, got a=" a ", b=" b))))

(defn -main [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args cli-options)
        {:keys [a b V0 mu steps Emin Emax out layers]} options]
    (when (:help options) (exit! 0 (usage summary)))
    (when (seq errors) (exit! 1 (str (str/join "\n" errors) "\n\n" (usage summary))))
    
    (let [energies (generate-energy-grid Emin Emax steps)]
      (if (some? layers)
        ;; Multilayer transfer-matrix mode
        (process-multilayer layers energies mu out)
        ;; Single KP closed-form mode
        (do
          (validate-single-kp-params a b)
          (process-single-kp a b V0 mu energies out))))))
