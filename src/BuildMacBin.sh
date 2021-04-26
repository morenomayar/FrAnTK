#replace zcat with gzcat in python scripts
for i in addBams.py getD.py getDNT.py getF4.py getF4NT.py getf3.py getf3NT.py getF4Ratio.py getF4RatioNT.py getF4subtr.py getF4subtrNT.py getPWdist.py getPWdistNT.py getDtrip.py getDtripNT.py getEnhD.py getEnhDNT.py getDstrat.py getDstratNT.py getDstrat2.py getDstrat2NT.py Freqs2Treemix.py CheckContigOrder.R
do
perl -pe 's/\[\"zcat\"/\[\"gzcat\"/g;' ../bin/$i > ../bin/tmp
mv ../bin/tmp ../bin/$i
done

cd ../bin
chmod 755    addBams.py getD.py getDNT.py getF4.py getF4NT.py getf3.py getf3NT.py getF4Ratio.py getF4RatioNT.py getF4subtr.py getF4subtrNT.py getPWdist.py getPWdistNT.py getDtrip.py getDtripNT.py getEnhD.py getEnhDNT.py getDstrat.py getDstratNT.py getDstrat2.py getDstrat2NT.py Freqs2Treemix.py CheckContigOrder.R
