#check if python exists
$pythonpath -c "for r in range(10): print(r)" > /dev/null 2>/dev/null
if [ $? -ne 0 ]
then
echo "$pythonpath"" is not a valid python path"
exit 1
fi

#check if perl exists
$perlpath -e 'for($i=0; $i<10; $i++){print $i;}' > /dev/null 2>/dev/null
if [ $? -ne 0 ]
then
echo "$perlpath"" is not a valid perl path"
exit 1
fi

#check if Rscript exists
$Rpath -e 'for(i in 1:10){print(i)}' > /dev/null 2>/dev/null
if [ $? -ne 0 ]
then
echo "$Rpath"" is not a valid R path"
exit 1
fi

#check if plink calls plink
a=`plink --version --noweb | perl -e '$l=<>; if($l =~ /PLINK v(\d.\d+)\S+/){$v=$1; if($v >= 1.9){print "0\n"; }}else{print "1\n"; }'`
if [ $a -ne 0 ]
then
echo "No plink. plinkv1.9 or higher should be available by just typing plink. Please set the \$PATH variable accordingly. "
exit 1
fi

#check if samtools calls samtools

a=`samtools --version | perl -e '$l=<>; if($l =~ /samtools (\d.\d+)/){$v=$1; if($v >= 1){print "0\n"; }}else{print "1\n"; }'`
if [ $a -ne 0 ]
then
echo "No samtools. samtools1 or higher should be available by just typing samtools. Please set the \$PATH variable accordingly. "
exit 1
fi

#create bin directory
mkdir ../bin

#build python scripts
for i in *.py
do
echo '#!'$pythonpath | cat - $i > ../bin/$i
done

#build perl scripts
for i in *.pl
do
echo '#!'$perlpath | cat - $i > ../bin/$i
done

#build R scripts
for i in *.R
do
echo '#!'$Rpath | cat - $i > ../bin/$i
done

#check if R libraries are installed
#doParallel
$Rpath -e 'library(doParallel)' > /dev/null 2>/dev/null
if [ $? -ne 0 ]
then
echo "doParallel R lib not installed for $Rpath. "
exit 1
fi

#ggplot2
$Rpath -e 'library(ggplot2)' > /dev/null 2>/dev/null
if [ $? -ne 0 ]
then
echo "ggplot2 R lib not installed for $Rpath. "
exit 1
fi

#gplots
#$Rpath -e 'library(gplots)' > /dev/null 2>/dev/null
#if [ $? -ne 0 ]
#then
#echo "gplots R lib not installed for $Rpath. "
#exit 1
#fi

#grid
$Rpath -e 'library(grid)' > /dev/null 2>/dev/null
if [ $? -ne 0 ]
then
echo "grid R lib not installed for $Rpath. "
exit 1
fi

cd ../bin
chmod 755 *

echo "Done. "
echo "IMPORTANT: make sure to add the new bin dir to your \$PATH variable. For example, you can add the following line to your ~/.bashrc"
echo "export PATH="$PWD":\$PATH"
