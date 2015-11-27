pro annular_PCA_by_lintrick, nod_beam, x_cen, y_cen, n_PCA, side

if side eq 0 then side='S'
if side eq 1 then side='D'

;////////////////////////////////////////////////////////
;Make polar and cartisian coordinate systems
;////////////////////////////////////////////////////////

X_DIM=300l
Y_DIM=300l

;polar
radius=fltarr(X_DIM,Y_DIM)
azimuth=fltarr(X_DIM,Y_DIM)

for i=0, X_DIM-1 do begin
   for j=0, Y_DIM-1 do begin
      ;radius        
      radius[i,j]=sqrt((i-x_cen)^2.+(j-y_cen)^2.)
      ;azimuth
      if i ne x_cen OR j ne y_cen then begin
         if i le x_cen AND j gt x_cen then azimuth[i,j]=atan((x_cen-i)/(j-y_cen))
         if i lt x_cen AND j le x_cen then azimuth[i,j]=atan((y_cen-j)/(x_cen-i))+!pi/2.
         if i ge x_cen AND j lt x_cen then azimuth[i,j]=atan((i-x_cen)/(y_cen-j))+!pi
         if i gt x_cen AND j ge x_cen then azimuth[i,j]=atan((j-y_cen)/(i-x_cen))+3.*!pi/2.
      endif
   endfor
endfor
azimuth*=360./(2*!pi) ;convert to degrees
stop

;cartesian
x_grid=findgen(X_DIM)#(fltarr(Y_DIM)+1.)
y_grid=reverse(findgen(Y_DIM))##(fltarr(X_DIM)+1)

;reform to 1-dimension
radius_1d=reform(radius,X_DIM*Y_DIM)
azimuth_1d=reform(azimuth,X_DIM*Y_DIM)
x_grid_1d=reform(x_grid,X_DIM*Y_DIM)
y_grid_1d=reform(y_grid,X_DIM*Y_DIM)

;////////////////////////////////////////////////////////
;Get the data
;////////////////////////////////////////////////////////

if nod_beam eq 0 and side eq 'D' then begin
spawn, 'ls ../../processed_data/sat/cxsxbcponlm/DX_0/*.fits', filenames
dir='../../processed_data/sat/cxsxbcponlm/DX_0/'
endif
if nod_beam eq 1 and side eq 'D' then begin
spawn, 'ls ../../processed_data/sat/cxsxbcponlm/DX_1/*.fits', filenames
dir='../../processed_data/sat/cxsxbcponlm/DX_1/'
endif
if nod_beam eq 0 and side eq 'S' then begin
spawn, 'ls ../../processed_data/sat/cxsxbcponlm/SX_0/*.fits', filenames
dir='../../processed_data/sat/cxsxbcponlm/SX_0/'
endif
if nod_beam eq 1 and side eq 'S' then begin
spawn, 'ls ../../processed_data/sat/cxsxbcponlm/SX_1/*.fits', filenames
dir='../../processed_data/sat/cxsxbcponlm/SX_1/'
endif

n_images=n_elements(filenames)
data=fltarr(X_DIM,Y_DIM,n_images)

;////////////////////////////////////////////////////////
;start a loop to do region-by-region
;////////////////////////////////////////////////////////
n_r=140

r_subt_in=findgen(n_r)*1.+3.
r_subt_out=findgen(n_r)*1.+3.+1.

r_opt_in=r_subt_in-7.
r_opt_out=r_subt_in+9.

r_opt_in[where(r_opt_in lt 6)]=5

;////////////////////////////////////////////////////////
;Do the whole thing in loops for memory allocation
;////////////////////////////////////////////////////////

combined_median=fltarr(X_DIM,Y_DIM)
combined_data=fltarr(X_DIM,Y_DIM,n_images)

for l=0, n_r-1 do begin
    print, l
       
    ;actually get the data here since we keep writing over the array
    for h=0, n_images-1 do data[*,*,h]=mrdfits(strcompress(filenames[h], /remove_all), /silent)
    for h=0, n_images-1 do data[*,*,h]-=median(data[*,*,h])

    ;////////////////////////////////////////////////////////
    ;construct a masked region
    ;////////////////////////////////////////////////////////
       
    opt_mask=fltarr(X_DIM*Y_DIM)
    opt_mask[where(radius_1d lt r_opt_in[l])]=1
    opt_mask[where(radius_1d ge r_opt_out[l])]=1

    n_pix_opt_masked=total(opt_mask)

    sub_mask=fltarr(X_DIM*Y_DIM)
    sub_mask[where(radius_1d lt r_subt_in[l])]=1
    sub_mask[where(radius_1d ge r_subt_out[l])]=1

    n_pix_sub_masked=total(sub_mask)
       
    ;////////////////////////////////////////////////////////
    ;make matrix of non-opt_masked pixels for PCA
    ;////////////////////////////////////////////////////////

    A=fltarr(X_DIM*Y_DIM-n_pix_opt_masked,n_images)

    for h=0, n_images-1 do begin
        data_1d=reform(data[*,*,h],X_DIM*Y_DIM)
        A[*,h]=data_1d[where(opt_mask eq 0)]
    endfor
       
    ;////////////////////////////////////////////////////////
    ;Find the principal components (eigenvectors)
    ;////////////////////////////////////////////////////////
       
    mean_image=total(A,2)/n_images
    A=A-rebin(mean_image,X_DIM*Y_DIM-n_pix_opt_masked,n_images) ;subtract off the mean image
       
    ;linear algebra desribed on the wiki page for eigenfaces
    AT_A=matrix_multiply(A,A, /atranspose)
    eigenvalues=eigenql(AT_A, eigenvectors=eigenvectors, /double)
    eigenvectors=A#eigenvectors

    ;normalize the eigenvectors
    for h=0, n_images-1 do eigenvectors[*,h]/=sqrt(total(eigenvectors[*,h]^2.))
       
    ;////////////////////////////////////////////////////////
    ;Fit the individual images and reform to 2d
    ;////////////////////////////////////////////////////////
       
    ;loop through image by image
    for h=0, n_images-1 do begin

        ;fit the image
        PCA_image_1d=fltarr(X_DIM*Y_DIM-n_pix_opt_masked)

        for k=0, N_PCA-1 do begin
            coefficient=total(eigenvectors[*,k]*A[*,h])   
            PCA_image_1d+=eigenvectors[*,k]*coefficient
        endfor

        ;subtract from original
        PCA_subtracted_image=A[*,h]-PCA_image_1d

        ;insert the opt_masked pixels back in
        PCA_subtracted_with_opt_mask_1d=fltarr(X_DIM*Y_DIM)
        PCA_subtracted_with_opt_mask_1d[where(opt_mask eq 0)]=PCA_subtracted_image

        ;only save the subtraction region (not the opt region)
        PCA_subtracted_image*=(1.-sub_mask)

        ;reform to 2-d
        data[*,*,h]=reform(PCA_subtracted_with_opt_mask_1d,X_DIM,Y_DIM)
    endfor
       
    ;////////////////////////////////////////////////////////
    ;rotate and combine
    ;////////////////////////////////////////////////////////

    restore, strcompress('../../processed_data/sat/cxsxbcponlm/'+side+'X_'+string(fix(nod_beam))+'_PA.sav', /remove_all)

    for h=0, n_images-1 do begin
        data[*,*,h]=rot(data[*,*,h],combined_images_PA[h], 1.0, x_cen, y_cen, /interp, /pivot, cubic=-0.5, missing=0.0)
    endfor

    ;median combine and save
    sub_mask_2d=reform(sub_mask, X_DIM,Y_DIM) 
    combined_median+=median(data,dimension=3)*(1-sub_mask_2d) ;only save the subtraction region, not the optimization region

    save_filename=strcompress('../../processed_data/sat/PCA/'+side+'X_'+string(fix(nod_beam))+'/median_annular_'+string(fix(n_PCA))+'.fits', /remove_all)
    mwrfits, combined_median, save_filename, /create

    ;   combined_data+=data
    ;ind_save_filename=strcompress('../../processed_data/sat/PCA/'+side+'X_'+string(fix(nod_beam))+'/individual_images/frame_'+string(h, format='(I03)')+'.fits', /remove_all)
    ;mwrfits, combined_data[*,*,h], ind_save_filename, /create

    endfor
end
