% New lensing code
% lensing.m
function lensing(s)
M = 1.99e30   ; % Mass of the Sun
c = 3e8       ; % Speed of light
G = 6.67e-11  ; % Newton's gravitational constant
Kp= 3.086e19  ; % A kiloParsec

Dol = 1 ; Dos = 2 ; % Distances in kpc

Dls = Dos - Dol ;

Sigma_cr = (c^2 / (4*pi*G) ) * Dos / ( Dol * Dls ) / Kp ; % Critical density

Er = sqrt( 4*G*M/c^2 * Dls * Dol / Dos * Kp ) ; % ER in lens plane

% Box size will be 8x8 Er

npix = 512 ;

den = zeros(npix,npix) ;
krn = zeros(npix,npix) ;

pix = 8 / npix ; % Pixel size in Ers
lpix= pix * Er ; % Pixel size in physical units (m)

den_pix = M / ( lpix * lpix ) ;

den(npix/4,npix/4) = den_pix / Sigma_cr ;

for i=1:npix
    for j=1:npix
        
        r = sqrt( (i-npix/2)^2 + (j-npix/2)^2 ) * pix ;
        
        if ( r > 0 )
            krn(i,j) = log( r ) ;
        else
            krn(i,j) = 0 ;
        end
    end
end

krn = circshift( krn , [-npix/2+1 -npix/2+1] );

den_fft = fft2( den ) ;
krn_fft = fft2( krn ) ;

den_fft = den_fft .* krn_fft ; 

den = real ( ifft2( den_fft ) ) ;

clear krn 
clear krn_fft
clear den_fft

td = den(1:npix/2,1:npix/2) * (pix * pix / pi ) ;

clear den

% put source at y = (s(1) , s(2))
for k=1:size(s,1)
x = zeros(npix/2) ;
y = zeros(npix/2) ;

for i = 1:npix/2
        
    for j=1:npix/2
        
        x(i,j) = (i-npix/4) * pix ;

        y(i,j) = (j-npix/4) * pix ;
        
        td(i,j) = (1/2)*( ( x(i,j)-s(k,1) )^2 + ( y(i,j)-s(k,2) )^2 ) - td(i,j) ;
        
    end 
end

figure('visible','on')
subplot(1,2,1)

contour(x,y,td,300)
hold on

xim1 = (hypot(s(k,1),s(k,2)) + sqrt( hypot(s(k,1),s(k,2))^2 + 4 ) )/2. ;
xim2 = (hypot(s(k,1),s(k,2)) - sqrt( hypot(s(k,1),s(k,2))^2 + 4 ) )/2. ;

plot([xim1 xim2],[0 0],'or')
axis square equal

subplot(1,2,2)
T = td(1:npix/2,npix/4) ;
xp = (-npix/4+1:npix/4)*pix ;
plot(xp,T)
hold on

% This should be the analytic result

t1 = (1/2)*( xim1 - hypot(s(k,1),s(k,2)) )^2 - log( abs(xim1) ) ;
t2 = (1/2)*( xim2 - hypot(s(k,1),s(k,2)) )^2 - log( abs(xim2) ) ;

plot([-npix/4*pix npix/4*pix],[t1 t1],'r')
plot([-npix/4*pix npix/4*pix],[t2 t2],'r')

end % for
end % function