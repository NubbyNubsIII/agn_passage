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


    cscres_table=pd.read_csv('/home/AGN_proj/agn_git/cscresults_nodup.tsv', sep='\t')
    cscobj_table=pd.read_csv('/home/AGN_proj/agn_git/csc_table.csv', sep=',')

  

    for i in range(len(par_coords_list)):

        current_par=par_names[i]

        try: 
            #load picture of par
            pic_hdulst=fits.open('/home/AGN_proj/agn_storage/big_pix/'+current_par+'_comb_drz_sci.fits')
            pic_header=pic_hdulst[0].header
            print(pic_header)
#            #get the wcs from the header

            jw_wcs=wcs.WCS(pic_header)
            data=pic_hdulst[0].data
            
            im_scale=np.abs(pic_header['CD1_1']) #degrees per pixel

        except FileNotFoundError:
            print("No Associated (Pic) Fits file for "+ current_par)
            continue


        try: 
            hdulst=fits.open('/home/AGN_proj/agn_storage/spec_fits/'+current_par+'_speccat.fits')
            tabhdu=hdulst[1]
            jw_table=tabhdu.data
        except FileNotFoundError:
            print("No Associated (table) Fits file for "+ current_par)
            continue

        print(current_par+'\n')

        cscobj_table_par=cscobj_table[cscobj_table['par_name'] == current_par]
        obj_coords_par=coord.SkyCoord(ra=cscobj_table['ra'], dec=cscobj_table['dec'], unit=(u.deg, u.deg))

        jw_coords=coord.SkyCoord(ra=jw_table['ra'],dec=jw_table['dec'],unit=(u.deg, u.deg))

        #now we want to loop throught the chandra objects in the par, and use the chan filenames to get the region image, and load that
        #then we 

        for source in range(len(cscobj_table_par)):

            source_info=cscobj_table_par[source]
            source_coords=obj_coords_par[source]

            #load the chandra region image

            chan_filename=source_info['chan_filename']

            try:
                chan_hdulst=fits.open('/home/AGN_proj/agn_storage/spec_fits/'+chan_filename)
                chan_header=chan_hdulst[0].header
                chan_data=chan_hdulst[0].data
                ch_wcs=wcs.WCS(chan_header)
            except FileNotFoundError:
                print("No Associated (chan) Fits file for "+ current_par + " " + source_info['name'])
                continue

            criteria=jw_coords.separation(source_coords) < source_info['error_ellipse_r0']*1.2*u.arcsec
            
            candidates= jw_table[criteria]
            cand_coords=jw_coords[criteria]
            #candidates is a table of all the candidates in the region, with their coordinates
            #cand_coords is a SkyCoord object of the candidates in the region
            if len(candidates)>0:
                print('\n Chandra Name: ' + source_info['name'])
                print('Significance: ' + str(source_info['significance']))
                print('Candidates:')
                print(candidates['id'])
#left to do:
                #make a plot of the chandra region image, with the candidates overlaid
                #and the error ellipse overlaid
                #and the PASSAGE JW image as either a rainbow contour plot (chand grayscale imshow)
                #or a different colored imshow (red and green, or blue) with opacity 0.5 for each

                #USE astropy.visualization.wcsaxes/overlay tutorials available instead
                #of what used for this last time.

                #basic idea - use WCS for the axes object - this is the chandra coords
                # then, use transform for the contour (JW plot)
                # maybe can also use transform for the double imshow version? 

                #maybe keep the original plot with just the ellipse and jw too for good context
                #would be awesome to have these side by side

                #dont worry about Cutout2D for JW image for now - but might be more performant?
                #plotting the chandra image as a
                x_jw,y_jw=w.world_to_pixel(cand_coords)
                id_jw=candidates['id']

                smaj_ax=source_info['error_ellipse_r0']
                smin_ax=source_info['error_ellipse_r1']
                ell_ang=source_info['error_ellipse_angle'] #angl of ellipse measured from north->east
                
                im_center=w.world_to_pixel(source_coords)
                im_edge_len=smaj_ax*3/(3600*im_scale)

                #plt.figure() JJUST CHECKING DISTRIBUTION FOR COLORSCALING
                #data_list=np.ravel(data)
                #plt.hist(data_list[np.logical_and(data_list>0,data_list<0.025)], bins=1000)
                #plt.show()

                fig=plt.figure(figsize=(8,8))

                 #make axes object
                ax=fig.add_axes((1.5/8,1.5/8,5/8,5/8),projection=w)

                overlay=ax.get_coords_overlay('icrs')
                overlay.grid(color='white',ls='dotted')
                overlay[0].set_axislabel('Right Ascension (degrees, J2000)') #label the overlay
                overlay[1].set_axislabel('Declination (degrees, J2000)')
                ax.set_xlabel("Right Ascension") #label the other sides for clarity
                ax.set_ylabel("Declination")
                ax.set_title(current_par+": CSC "+source_info['_2CXO']+", Sig: "+str(source_info['S_N']),
                             y=0.95, bbox=dict(boxstyle="square,pad=0.3",
                                               fc='lightblue'))
                
                ell=Ellipse(im_center,2*smin_ax/(3600*im_scale),
                            2*smaj_ax/(3600*im_scale),angle=-ell_ang,
                            fill=False,color='red')
                ax.add_patch(ell)

                ax.set_xlim(im_center[0]-im_edge_len/2,im_center[0]+im_edge_len/2)
                ax.set_ylim(im_center[1]-im_edge_len/2,im_center[1]+im_edge_len/2)

                ax.scatter(im_center[0],im_center[1],marker='+',color='red')

                

                for i, txt in enumerate(id_jw):
                    ax.annotate(txt,(x_jw[i],y_jw[i]),xytext=(x_jw[i]+10,y_jw[i]),
                                bbox=dict(boxstyle="square,pad=0.1",
                      fc="white", ec="black"))
                ax.scatter(x_jw,y_jw,marker='*',s=10,color='cyan')    


                cnorm = mpl.colors.SymLogNorm(0.0005,0.1)

                cmap = 'gray' #use magma color map
                ### show the data on an image object, normalized by lognorm object
                im = ax.imshow(data, origin='lower', norm=cnorm, cmap=cmap)
                
                chan_for_file=source_info['_2CXO'].replace("+","_")
                plt.savefig(current_par+"_CSC_"+chan_for_file+".png")
                


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

