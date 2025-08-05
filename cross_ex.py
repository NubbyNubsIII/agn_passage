import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.coordinates as coord
import numpy as np
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.wcs import WCS
import astropy.units as u
from astropy import wcs
from matplotlib.patches import Ellipse
import pandas as pd
from astropy.nddata import Cutout2D
import astropy.nddata.utils as utils

#STUFF TO DO 
#1 load big picture, get wcs from header
#2 convert all sky coords to wcs world to pixel
#3 do overlay stuff to get RA Dec plot nicely
#4 look at the inset plots code to plot just the inset, with diameter 3*r0
#5 figure out how to plot the ellipse (especially the rotation, need to add how far away declination is from north)


#functions:

#f

if __name__ == "__main__":
    
    # Place a big circle at each par TARGET coords
    # should be 3 arcminutes radius
    #get all chandra objects in that radius


    par_coords_import=np.loadtxt('./par_coords.csv',dtype='str',delimiter=',')
    par_coords_list=coord.SkyCoord(ra=par_coords_import[:,1],dec=par_coords_import[:,2],unit=(u.deg, u.deg))
    par_names=par_coords_import[:,0]


    cscres_table=pd.read_csv('/home/jacoblevine7/AGN_proj/agn_git/cscresults_nodup.tsv', sep='\t')
    cscobj_table=pd.read_csv('/home/jacoblevine7/AGN_proj/agn_git/csc_table.csv', sep=',')

  

    for i in range(len(par_coords_list)):

        current_par=par_names[i]

        try: 
            #load picture of par
            pic_hdulst=fits.open('/home/jacoblevine7/AGN_proj/agn_storage/big_pix/'+current_par+'_comb_drz_sci.fits')
            pic_header=pic_hdulst[0].header
#            #get the wcs from the header

            jw_wcs=wcs.WCS(pic_header)
            data=pic_hdulst[0].data
            
            im_scale=np.abs(pic_header['CD1_1']) #degrees per pixel

        except FileNotFoundError:
            print("No Associated (Pic) Fits file for "+ current_par)
            continue


        try: 
            hdulst=fits.open('/home/jacoblevine7/AGN_proj/agn_storage/spec_fits/'+current_par+'_speccat.fits')
            tabhdu=hdulst[1]
            jw_table=tabhdu.data
        except FileNotFoundError:
            print("No Associated (table) Fits file for "+ current_par)
            continue

        print(current_par+'\n')

        cscobj_table_par=cscobj_table[cscobj_table['par_name'] == current_par]
        print(f'Number of Chandra Objects in {current_par}: {len(cscobj_table_par)}')
        obj_coords_par=coord.SkyCoord(ra=cscobj_table_par['ra'], dec=cscobj_table_par['dec'], unit=(u.deg, u.deg))

        jw_coords=coord.SkyCoord(ra=jw_table['ra'],dec=jw_table['dec'],unit=(u.deg, u.deg))

        #now we want to loop throught the chandra objects in the par, and use the chan filenames to get the region image, and load that
        #then we 

        for source in range(len(cscobj_table_par)):

            source_info=cscobj_table_par.iloc[source]
            source_coords=obj_coords_par[source]

            #load the chandra region image

            chan_filename=source_info['chan_filename']
                        
            smaj_ax=source_info['error_ellipse_r0']
            smin_ax=source_info['error_ellipse_r1']
            ell_ang=source_info['error_ellipse_angle']


            try:
                chan_hdulst=fits.open('/home/jacoblevine7/AGN_proj/agn_storage/chan_pix/chan_pix_good/'+chan_filename)
                chan_header=chan_hdulst[0].header
                chan_data=chan_hdulst[0].data
                chan_wcs=wcs.WCS(chan_header)
            except FileNotFoundError:
                print("No Associated (chan) Fits file for "+ current_par + " " + source_info['name'])
                continue

            criteria=jw_coords.separation(source_coords) < np.max([source_info['error_ellipse_r0'].item()*1.3,2.7])*u.arcsec
            
            candidates= jw_table[criteria]
            cand_coords=jw_coords[criteria]
            #candidates is a table of all the candidates in the region, with their coordinates
            #cand_coords is a SkyCoord object of the candidates in the region
            if len(candidates)>-1:
                print('\n Chandra Name: ' + source_info['name'])
                print('Significance: ' + str(source_info['significance']))
                print('Candidates:')
                print(candidates['id'])
#left to do:
                #make shit look nice (standardize plot sizes, etc)
                #cull the bad ones
                #explore colorbar on top of right plot?
                #consider plotting the whole par with all the chandra objects too, to see which ones are actually in the field
                #need to figure out if the objects that appear not to be labelled are actually labelled or not?

                
                cnorm_jw=mpl.colors.SymLogNorm(0.0005,0.1)
                
                size=np.max([source_info['error_ellipse_r0'].item()*3,6])*u.arcsec # put a minimum size of 6 arcsec

                try:
                    cutout=Cutout2D(data,position=source_coords,size=size,wcs=jw_wcs)
                except utils.NoOverlapError:
                    print(f'Cutout for {source_info["name"]} is out of bounds of the image.')
                    continue
                im_center=source_coords.to_pixel(cutout.wcs)
                x_jw,y_jw=cand_coords.to_pixel(cutout.wcs)
                id_jw=candidates['id']
                #make a colormap to use for the overlay data
                base_cmap = plt.cm.hot
                new_cmap = mpl.colors.ListedColormap(list(base_cmap(np.linspace(0, 1, 256))))
                new_cmap.set_under(color=(1,0,0,0))  # Set the color for out_of_bounds values to transparent red

                fig=plt.figure(figsize=(20,5))

                gs=fig.add_gridspec(1,4,hspace=0,wspace=0)
                ax1 = fig.add_subplot(gs[0, 0], projection=cutout.wcs)
                ax2 = fig.add_subplot(gs[0, 1], projection=cutout.wcs)  
                ax3 = fig.add_subplot(gs[0, 2], projection=chan_wcs)   

                ax1.tick_params(axis='x',top=False, bottom=True, labelbottom=True)
                ax1.tick_params(axis='y',right=False, left=True,labelleft=True)
                ax1.set(xlabel='RA', ylabel='Dec')
                ax1.set_title('Combined Image')

                ax2.tick_params(axis='x',top=False, left=False, right=False,labelleft=False,labelbottom=False)
                ax2.tick_params(axis='y',top=False, left=False, right=False,labelleft=False,labelbottom=False)
                ax2.set(xlabel='RA', ylabel='Dec')
                ax2.set_title('JWST Field Cutout')

                ax3.tick_params(axis='x',top=False, labelbottom=True, bottom=True)
                ax3.tick_params(axis='y', right=False, labelleft =False, labelright=True,left=False)
                ax3.set(xlabel='RA',ylabel=' ')
                ax3.set_title('Chandra Region Image')


                ax1.set_xlim(0, cutout.data.shape[1])
                ax1.set_ylim(0, cutout.data.shape[0])


                ax1.imshow(cutout.data,origin='lower',cmap='gray',norm=cnorm_jw)
                ax1.imshow(chan_data,origin='lower', cmap=new_cmap, vmin=0,transform=ax1.get_transform(chan_wcs), alpha=0.3)
                ax1.grid(color='white', alpha=0.2, linestyle='--')

                ax2.imshow(cutout.data,origin='lower',cmap='gray',norm=cnorm_jw)
                ax2.text(cutout.data.shape[0]/2, -cutout.data.shape[1]/20, s=str(source_info['name']), fontsize=12, color='black', ha='center', va='center')
                ax2.text(cutout.data.shape[0]/2, -cutout.data.shape[1]/8, s=str(source_info['par_name']), fontsize=12, color='black', ha='center', va='center')

                ell=Ellipse(im_center,2*smin_ax/(3600*im_scale),
                                            2*smaj_ax/(3600*im_scale),angle=-ell_ang,
                                            fill=False,color='red')
                ax2.add_patch(ell)
                ax2.scatter(im_center[0],im_center[1],marker='+',color='red',s=20)

                for i, txt in enumerate(id_jw):
                    ax2.annotate(txt,(x_jw[i],y_jw[i]),xytext=(x_jw[i]+5,y_jw[i]),
                                bbox=dict(boxstyle="square,pad=0.1",
                        fc="white", ec="black"))
                    ax2.scatter(x_jw,y_jw,marker='*',s=8,color='cyan')


                ax3.imshow(chan_data,origin='lower', cmap='hot')
                ax3.grid(color='white', alpha=0.2, linestyle='--')
                ax3.text(chan_data.shape[0]-chan_data.shape[0]/6, chan_data.shape[1]-chan_data.shape[1]/15, s='Sig: '+str(np.round(source_info['significance'],3)),
                    fontsize=10, ha='center', va='center',bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="black"))
                #ax3.annotate('Sig: '+str(np.round(source['significance'].item(),4)),(chan_data.shape[0]-chan_data.shape[0]/6, chan_data.shape[1]-chan_data.shape[1]/6),xytext=(x_jw[i]+5,y_jw[i]),
                #              bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="black"))


                
                fig.savefig('/home/jacoblevine7/AGN_proj/agn_storage/chan_results/'+current_par+"_CSC_"+source_info['name'].split(' ',1)[1]+".png")
                


            else:
                continue
                print('Chandra Name: ' + source_info['_2CXO'])
                print('No Candidates')
                print('Closest Object vs edge of confident chandra area')
                print(np.min(par_coords_list.separation(chan_coords[source])))
                print('is greater than')
                print(source_info['r0']*u.arcsec)


            
            #use a boolean vectorized thingy to make a mask for the skycoords outside of error
            #i think that skycoord has a separation method you can use
            #then if you apply the mask, then maybe it will just make the correct selection of
            #par table. Still need to figure out how to get the chandra information into metadata
            # aka in the header. should be easy to do.

        #we want to somewhere have a method for ordering by significance
            #and prioritize ones with less potential objects?
    

    #now we have all the chandra information for the obj in the par range
    #now we wanna load the speccat.fits file for that par

