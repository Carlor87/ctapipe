"""
Simple pipeline using ctapipe to read HESSfits files generated after image cleaning
using HAP

Author: Carlo Romoli - MPIK

"""
get_ipython().magic(u'pylab')
from ctapipe.io import event_source
from ctapipe.instrument import CameraGeometry
from ctapipe.image import hillas, tailcuts_clean
from ctapipe.reco import hillas_intersection  ## Dan's module
from ctapipe.coordinates import *
from ctapipe.reco import HillasReconstructor
from ctapipe.visualization import CameraDisplay
import astropy.units as u
import matplotlib.pyplot as plt
#import numpy as np
from tqdm import tqdm
from astropy.table import Table

from IPython import display

# load immediately the canvas with the 4 subplots
# fig, axs = plt.subplots(1, 4, figsize=(15, 4), sharey=True, sharex=True)

def draw_several_cams_hillas(event,hillas):
    """
    Function to plot camera images on a canvas.
    :param event: event from the generator
    :param hillas: hillas parameters for the overlay
    :return: nothing
    """
    cmaps = 'jet'
    ncams = len(event.inst.subarray.tel)
    fig, axs = plt.subplots(1, 4, figsize=(15, 4), sharey=True, sharex=True)
    for ii in range(ncams):
        disp = CameraDisplay(
            event.inst.subarray.tel[ii+1].camera,
            ax=axs[ii],
            title="CT{}".format(ii + 1),
        )
        disp.cmap = cmaps
        try:
           disp.image = event.dl1.tel[ii+1].image
        except KeyError:
            continue
        if hillas[ii+1]:
            disp.overlay_moments(
                hillas[ii+1], color='pink', lw=3, with_label=False
            )
        else:
            continue


def draw_several_cams(event):
    """
    Function to plot camera images on a canvas.
    :param event: event from the generator
    :return: nothing
    """
    cmaps = 'jet'
    ncams = len(event.inst.subarray.tel)
    fig, axs = plt.subplots(1, 4, figsize=(15, 4), sharey=True, sharex=True)
    for ii in range(ncams):
        disp = CameraDisplay(
            event.inst.subarray.tel[ii+1].camera,
            ax=axs[ii],
            title="CT{}".format(ii + 1),
        )
        disp.cmap = cmaps
        try:
           disp.image = event.dl1.tel[ii+1].image
        except KeyError:
            continue


def PlotNewEv():
    event = next(gsource)
    draw_several_cams(event)
    plt.show()


def PlotSummedImages(event,mask,hillas):
    # plt.figure()
    # geom = event.inst.subarray.tel[-1].camera
    disp = CameraDisplay(geom,title="HESS I")
    disp.cmap = 'jet'
    image=0
    for ii in event.dl0.tels_with_data:
        image += event.dl1.tel[ii].image*mask[ii]
        disp.image = image
        disp.overlay_moments(hillas[ii], color='pink', lw=3, with_label=False,keep_old=True)


def PlotSummedImages_fromList(ev,image,mask,hillas):
    #plt.figure()
    #geom = event.inst.subarray.tel[-1].camera
    plt.cla()
    disp = CameraDisplay(geom,title="HESS I",ax=ax)
    disp.cmap = 'jet'
    sum_image=0
    for ii in image[ev].keys():
        sum_image += image[ev][ii]*mask[ev][ii]
        disp.image = sum_image
        disp.overlay_moments(hillas[ev][ii], color='pink', lw=3, with_label=False,keep_old=True)


def PlotSeveralCams_fromList(ev,image,hillas):
    cmaps = 'jet'
    ncams = 4
    fig, axs = plt.subplots(1, 4, figsize=(15, 4), sharey=True, sharex=True)
    plt.cla()
    for ii in range(ncams):
        disp = CameraDisplay(
            geom,
            ax=axs[ii],
            title="CT{}".format(ii + 1),
        )
        disp.cmap = cmaps
        try:
            disp.image = image[ev][ii+1]
        except KeyError:
            continue
        if hillas[ev][ii+1]:
            disp.overlay_moments(hillas[ev][ii + 1], color='pink', lw=3, with_label=False)
        else:
            continue

def read_pointing(filename):
    '''
    Function to read the pointing correction file
    The pointing correction file is an ascii table
    compiled according to the ecsv format.
    :return: pointing table
    '''
    table = Table(filename,format='ascii.ecsv')



if __name__ == "__main__":

    camdef=Table.read("/home/cromoli/Documents/hessfits/pks2155flarerun_wImages_test.hap/chercam.fits.gz",format='fits') #read the table with the pixel position
    # separate the different cameras in different tables
    camdef1=camdef[:960]
    camdef2=camdef[960:1920]
    camdef3=camdef[1920:2880]
    camdef4=camdef[2880:3840]
    geom = CameraGeometry('HESSI',camdef3['PIX_ID'],camdef3['PIX_POSX']*u.m,camdef3['PIX_POSY']*u.m,camdef3['PIX_AREA']*u.m*u.m,pix_type='hexagonal')

    ShowerList=[]
    ImageList=[]
    MaskList=[]
    HillasList=[]

    count =0

    with event_source("/home/cromoli/Documents/hessfits/pks2155flarerun_wImages_test.hap/run_00033746_televentlist_img.fits") as source:
        gsource = (x for x in source)
        # hillasReco = HillasReconstructor() ## method 1
        hillasReco = hillas_intersection.HillasIntersection()  ## method 2
        for i in tqdm(range(100)):
#       for event in source:       # to read the entire run
            event = next(gsource)   # load a single event as a test
            HillasIm=dict()
            NomHillasIm=dict()
            ExtendedImage=dict()   # only for presentation purposes
            mask0510=dict()
            telx=dict()
            tely=dict()
            telphi=dict()
            teltheta=dict()
            ''' !!!!this loop is slow!!!! '''
            for j in event.dl0.tels_with_data:
                foc=event.inst.subarray.tel[j].optics.equivalent_focal_length
                conv = u.rad/foc
                # test to check the reconstruction module
#                telphi[j]=event.pointing.azimuth*u.deg
#                teltheta[j]=event.pointing.altitude*u.deg
                ExtendedImage[j]=event.dl1.tel[j].image
                mask0510[j] = tailcuts_clean(event.inst.subarray.tel[j].camera,event.dl1.tel[j].image,10,5) #implement a (5,10) cleaning
                HillasIm[j]=hillas.hillas_parameters(event.inst.subarray.tel[j].camera,event.dl1.tel[j].image * mask0510[j])
                NomHillasIm[j]=HillasIm[j]._replace(cen_x=HillasIm[j].cen_x*conv,
                                                cen_y=HillasIm[j].cen_x*conv,
                                                length=HillasIm[j].length*conv,
                                                width=HillasIm[j].width*conv,
                                                r=HillasIm[j].r*conv)
                telx[j]=event.inst.subarray.pos_x[j - 1]
                tely[j]=event.inst.subarray.pos_y[j - 1]

            ImageList.append(ExtendedImage)  ## only for presentation purposes
            MaskList.append(mask0510)
            HillasList.append(HillasIm)
            ''' reconstructor in the nominal system of the cameras '''
            array_pointing = HorizonFrame(alt=event.pointing.altitude*u.deg, az=event.pointing.azimuth*u.deg)
            shower = hillasReco.predict(NomHillasIm, telx, tely, array_pointing)
            # shower = hillasReco.predict(hillas_dict=HillasIm,inst=event.inst,tel_phi=telphi,tel_theta=teltheta)
            ShowerList.append(shower)
        print("End process!")

    fig=plt.figure()
    ax=plt.gca()
    PlotSummedImages_fromList(ev=10,image=ImageList,mask=MaskList,hillas=HillasList)
    PlotSeveralCams_fromList(ev=10,image=ImageList,hillas=HillasList)

    ## check some reconstructed features
    corx=[x['core_x']/u.m for x in ShowerList]
    cory=[x['core_y']/u.m for x in ShowerList]
    alt=[x['alt'][0]/u.deg for x in ShowerList]
    az=[x['az'][0]/u.deg for x in ShowerList]

    plt.figure()
    plt.hist2d(az,alt,bins=[np.linspace(155,210,50),np.linspace(75,90,50)])
    plt.xlabel("Azimuth angle [deg]")
    plt.ylabel("Altitude angle [deg]")

    plt.figure()
    plt.hist2d(corx,cory,bins=[np.linspace(-200,200,50),np.linspace(-200,200,50)])
    plt.show()

