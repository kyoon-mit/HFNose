for i in {50..300..50}
do cmsRun TICLAnalyzer.py inpath=electron_2026D60_HFNose outpath=electron_2026D60_HFNose process=electron pid=11 E=${i}
cmsRun TICLAnalyzer.py inpath=photon_2026D60_HFNose outpath=photon_2026D60_HFNose process=photon pid=22 E=${i}
cmsRun TICLAnalyzer.py inpath=charged_pion_2026D60_HFNose outpath=charged_pion_2026D60_HFNose process=charged_pion pid=211 E=${i}
done
