
## for max-min driving force and enzyme protein cost analysis 
__WLP__
```
python path\to\PathParser\Scripts\main1.py -o path\to\output -r path\to\WLP.tsv -eb ATP,ADP,Pi,NADH,NAD,NADPH,NADP -b 0.001,10 -w 12
```
__Bicycle_CO2__
```
python path\to\PathParser\Scripts\main1.py -o path\to\output -r path\to\Bicycle_CO2.tsv -eb ATP,ADP,Pi,NADH,NAD,NADPH,NADP,rFdx,oFdx,CO2 -b 0.001,10 -w 12 -a pfor:1
```
__Bicycle_formate__
```
python path\to\PathParser\Scripts\main1.py -o path\to\output -r path\to\Bicycle_formate.tsv -eb ATP,ADP,Pi,NADH,NAD,NADPH,NADP,Fm -b 0.001,10 -w 12 -a pfl:1
```
## for robustness and flux fold change analysis
__WLP__
```
python path\to\PathParser\Scripts\main2.py -o path\to\output -r path\to\WLP.tsv -eb ATP,ADP,Pi,NADH,NAD,NADPH,NADP -eo  ATP,ADP,Pi,NADH,NAD,NADPH,NADP -b 0.1,10 -f Ac -n 100 -w 123 -d no -p 30 -t no
```
__Bicycle_CO2__
```
python path\to\PathParser\Scripts\main2.py -o path\to\output -r path\to\Bicycle_CO2.tsv -eb ATP,ADP,Pi,NADH,NAD,NADPH,NADP,rFdx,oFdx -eo  ATP,ADP,Pi,NADH,NAD,NADPH,NADP,rFdx,oFdx -b 0.1,10 -i CO2 -n 100 -w 123 -d no -p 30 -t no
```
__Bicycle_formate__
```
python path\to\PathParser\Scripts\main2.py -o path\to\output -r path\to\Bicycle_formate.tsv -eb ATP,ADP,Pi,NADH,NAD,NADPH,NADP -eo  ATP,ADP,Pi,NADH,NAD,NADPH,NADP -b 0.1,10 -i Fm -n 100 -w 123 -d no -p 30 -t no
```
