
#!/bin/sh -vx

while read line
do
cp "$line" /home/jeffrey/Desktop/mass_8_sersic_2/
done < /home/jeffrey/Documents/dataproducts/mass_8_sersic_2.csv
