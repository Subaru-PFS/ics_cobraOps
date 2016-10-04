function output = getCentersRails(numOfRails)

% Create Rail Centers
centers1 = complex([0:8:(28*8)],0);
centers2 = centers1(1:end-1) + 8 * exp(1i * pi / 3);
centers = [centers1 centers2].' + 8 * exp(i*2*pi/3);
rails = centers;
for jj = 1:(numOfRails-1)
    rails = [rails; centers + jj * 16 * exp(1i * 2 * pi/3)];
end
output = rails;
end