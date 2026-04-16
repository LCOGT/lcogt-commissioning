rm *.png
mkdir plots
readmode="full_frame_lownoise"
sortby="filterlevel"
noisegainmef --readmode $readmode --sortby $sortby --makepng --minx  500 --maxx 1500 --miny 500 --maxy 1500  $@
for filename in *.png; do mv "$filename" "plots/${filename%%.*}_ll.png"; done;

noisegainmef --readmode $readmode --sortby $sortby --makepng --minx 2500 --maxx 3500 --miny 500 --maxy 1500  $@
for filename in *.png; do mv "$filename" "plots/${filename%%.*}_lr.png"; done;

noisegainmef --readmode $readmode --sortby $sortby --makepng --minx  500 --maxx 1500 --miny 2500 --maxy 3500  $@
for filename in *.png; do mv "$filename" "plots/${filename%%.*}_ul.png"; done;

noisegainmef --readmode $readmode --sortby $sortby --makepng --minx 2500 --maxx 3500 --miny 2500 --maxy 3500  $@
for filename in *.png; do mv "$filename" "plots/${filename%%.*}_ur.png"; done;
