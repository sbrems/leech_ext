pro make_skimage_coeffs
Kx = [[-2.1478925, 0.0058138110,  -7.6396687e-06,   2.5359596e-09], [1.0109149,  -2.3826537e-05,   2.8458629e-08,  -9.3206482e-12], [-2.1164521e-05, 5.3115381e-08,  -6.6315643e-11,   2.2888432e-14], [1.2983972e-08,  -4.1253977e-11,   5.1637044e-14,  -1.5988376e-17]]
outkx = fltarr(28)
outkx[0]=Kx[0,0]
outkx[1]=Kx[0,1]
outkx[2]=Kx[1,0]
outkx[3]=Kx[0,2]
outkx[4]=Kx[1,1]
outkx[5]=Kx[2,0]
outkx[6]=Kx[0,3]
outkx[7]=Kx[1,2]
outkx[8]=Kx[2,1]
outkx[9]=Kx[3,0]
outkx[11]=Kx[1,3]
outkx[12]=Kx[2,2]
outkx[13]=Kx[3,1]
outkx[17]=Kx[2,3]
outkx[18]=Kx[3,2]
outkx[24]=Kx[3,3]
print,'a'
print, outkx
Ky = [[9.2717864, 0.98776733,   4.3514612e-06,   9.3450739e-09], [-0.013617797,  -3.9526096e-05,   8.1204222e-08,  -5.2048768e-11], [1.1313247e-05,   6.7127301e-08,  -1.6531988e-10,   1.0656544e-13], [1.6283111e-09,  -2.7723216e-11,   8.2118035e-14,  -5.3695050e-17]]
outky = fltarr(28)
outky[0]=Ky[0,0]
outky[1]=Ky[0,1]
outky[2]=Ky[1,0]
outky[3]=Ky[0,2]
outky[4]=Ky[1,1]
outky[5]=Ky[2,0]
outky[6]=Ky[0,3]
outky[7]=Ky[1,2]
outky[8]=Ky[2,1]
outky[9]=Ky[3,0]
outky[11]=Ky[1,3]
outky[12]=Ky[2,2]
outky[13]=Ky[3,1]
outky[17]=Ky[2,3]
outky[18]=Ky[3,2]
outky[24]=Ky[3,3]
print, 'b'
print, outky
end
