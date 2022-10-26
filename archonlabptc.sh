rm *.png
mkdir plots

noisegainmef *.fits --readmode None --makepng --minx  500 --maxx 1500 --miny 500 --maxy 1500
for filename in *.png; do mv "$filename" "plots/a_${filename}"; done;

noisegainmef *.fits --readmode None --makepng --minx 2500 --maxx 3500 --miny 500 --maxy 1500
for filename in *.png; do mv "$filename" "plots/b_${filename}"; done;

noisegainmef *.fits --readmode None --makepng --minx  500 --maxx 1500 --miny 2500 --maxy 3500
for filename in *.png; do mv "$filename" "plots/c_${filename}"; done;

noisegainmef *.fits --readmode None --makepng --minx 2500 --maxx 3500 --miny 2500 --maxy 3500
for filename in *.png; do mv "$filename" "plots/d_${filename}"; done;
