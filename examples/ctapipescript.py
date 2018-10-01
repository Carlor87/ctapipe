"""
Simple pipeline using ctapipe to read HESSfits files generated after image cleaning
using HAP

Best usable with ipython

Author: Carlo Romoli - MPIK

"""
get_ipython().magic(u'pylab')
import astropy.units as u
import matplotlib.pyplot as plt
# copy module needed to have exact and independent copy of the original object
import copy
from ctapipe.io import event_source
from ctapipe.instrument import CameraGeometry
from ctapipe.image import hillas, tailcuts_clean
from ctapipe.reco import hillas_intersection  # Dan's module
from ctapipe.coordinates import *
from ctapipe.reco import HillasReconstructor  # alternative module for Hillas reconstruction
from ctapipe.visualization import CameraDisplay
# import numpy as np
from tqdm import tqdm
from astropy.table import Table
from ctapipe.calib import CameraCalibrator

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
    return table

def CT5CameraImage(ev,image,hillas):
    geom = CameraGeometry.from_name('HESS-II')
    cmaps = 'jet'
    plt.figure(figsize=(8,8))
    axs = plt.gca()
    disp = CameraDisplay(
            geom,
            ax=axs,
            title="CT{}".format(5),
        )
    disp.cmap = cmaps
    try:
        disp.image = image[ev][5]
    except KeyError:
        print('Dict error')
    if hillas[ev][5]:
        disp.overlay_moments(hillas[ev][5], color='pink', lw=3, with_label=False)
    else:
        print('no Hillas')


if __name__ == "__main__":

    camdef=Table.read("/home/cromoli/Documents/hessfits/pks2155flarerun_wImages_test.hap/chercam.fits.gz",format='fits') #read the table with the pixel position
    # separate the different cameras in different tables
    camdef1=camdef[:960]
    camdef2=camdef[960:1920]
    camdef3=camdef[1920:2880]
    camdef4=camdef[2880:3840]
    #geom = CameraGeometry(
    #    'HESSI',
    #    camdef3['PIX_ID'],
    #    np.array(camdef3['PIX_POSX']) * u.m,
    #    np.array(camdef3['PIX_POSY']) * u.m,
    #    np.array(camdef3['PIX_AREA']) * u.m * u.m,
    #    pix_type='hexagonal')  # Real camera
    geom = CameraGeometry.from_name('HESS-I') # only for MC

    ShowerList = []
    ImageList = []
    MaskList = []
    HillasList = []
    MCList = []

    count =0

    cal = CameraCalibrator(None, None, r1_product='HESSIOR1Calibrator', extractor_product ='FullIntegrator')

    # Real data source
    # with event_source(
    #       "/home/cromoli/Documents/hessfits/pks2155flarerun_wImages_test.hap/run_00033746_televentlist_img.fits.gz") as source:

    # MC source
    with event_source(
            "/home/cromoli/Documents/hessfits/gamma_0deg_0deg_run47444___phase2b2_desert-ws0_off0.5.simhess.gz") as source:
        # source.allowed_tels = [5]  # set allowed telescopes
        # gsource = (x for x in source)
        # hillasReco = HillasReconstructor() ## method 1
        hillasReco = hillas_intersection.HillasIntersection()  ## method 2
        # for i in tqdm(range(100)):
        for event in source:          # to read the entire run
            # event = next(gsource)   # load a single event as a test

            # dictionaries for storage purposes
            HillasIm = dict()
            NomHillasIm = dict()
            ExtendedImage = dict()   # only for presentation purposes
            mask0510 = dict()
            telx = dict()
            tely = dict()
            telphi = dict()
            teltheta = dict()
            mcdic = dict()

            cal.calibrate(event) # simple calibration for MC data

            for j in event.dl0.tels_with_data:
                foc = event.inst.subarray.tel[j].optics.equivalent_focal_length
                conv = u.rad/foc
                # test to check the reconstruction method 1
                # telphi[j]=event.mcheader.run_array_direction[0]*u.deg
                # teltheta[j]=event.mcheader.run_array_direction[1]*u.deg
                event.dl1.tel[j].image = event.dl1.tel[j].image[0]
                ExtendedImage[j] = event.dl1.tel[j].image
                mask0510[j] = tailcuts_clean(event.inst.subarray.tel[j].camera,
                                             event.dl1.tel[j].image,
                                             10,5)  # implement a (5,10) cleaning
                try:
                    HillasIm[j] = hillas.hillas_parameters(event.inst.subarray.tel[j].camera,
                                                           event.dl1.tel[j].image * mask0510[j])
                except hillas.HillasParameterizationError:
                    continue
                NomHillasIm[j] = copy.deepcopy(HillasIm[j])
                NomHillasIm[j].update(x=NomHillasIm[j].x*conv,
                                      y=NomHillasIm[j].y*conv,
                                      length=NomHillasIm[j].length*conv,
                                      width=NomHillasIm[j].width*conv,
                                      r=NomHillasIm[j].r*conv)
                telx[j] = event.inst.subarray.tel_coords.x[j - 1]
                tely[j] = event.inst.subarray.tel_coords.y[j - 1]

            ImageList.append(ExtendedImage)  # only for presentation purposes
            MaskList.append(mask0510)
            HillasList.append(HillasIm)
            # if len(event.dl0.tels_with_data) < 2:
            if len(NomHillasIm) < 2:
                continue
            ''' reconstructor in the nominal system of the cameras '''
            # array_pointing = HorizonFrame(alt=event.pointing.altitude*u.deg, az=event.pointing.azimuth*u.deg)
            try:
                array_pointing = HorizonFrame(alt=event.mcheader.run_array_direction[1].to('deg'),
                                              az=event.mcheader.run_array_direction[0].to('deg'))
                shower = hillasReco.predict(NomHillasIm, telx, tely, array_pointing)
                # shower = hillasReco.predict(hillas_dict=HillasIm,
                #                             inst=event.inst,
                #                             pointing_az=telphi,
                #                             pointing_alt=teltheta) # for method 1
            except ZeroDivisionError:
                continue
            mcdic = copy.deepcopy(event.mc)
            MCList.append(mcdic)
            ShowerList.append(shower)

        print("End process!")

    # PlotSeveralCams_fromList(94,ImageList,HillasList)
    # CT5CameraImage(94,ImageList,HillasList)

    # check some reconstructed features
    corx = [x['core_x'].value for x in ShowerList]
    cory = [x['core_y'].value for x in ShowerList]
    alt = [x['alt'][0].value for x in ShowerList]
    az = [x['az'][0].value for x in ShowerList]
    hmax = [x['h_max'].value for x in ShowerList]

    alt_mc = [x['alt'].to('deg').value for x in MCList]
    az_mc = [x['az'].to('deg').value for x in MCList]
    corx_mc = [x['core_x'].value for x in MCList]
    cory_mc = [x['core_y'].value for x in MCList]
    hmax_mc = [x['x_max'].value for x in MCList]
    ene_mc = [x['energy'].value for x in MCList]


    # plot core distance
    plt.figure()
    plt.hist2d(corx,cory,bins=[np.linspace(-200,200,50),np.linspace(-200,200,50)])
    plt.xlabel("Reconstructed Core X [m]")
    plt.ylabel("Reconstructed Core Y [m]")
    plt.show()

    # plot altitude comparison
    plt.figure()
    plt.plot(ene_mc, (-np.array(alt) + np.array(alt_mc))/np.array(alt_mc), 'o')
    plt.xlabel("True Energy [TeV]")
    plt.ylabel("(MC Altitude - Altitude)/(MC altitude)")
    plt.grid()
    plt.semilogx()

    # plot hmax comparison
    plt.figure()
    plt.plot(ene_mc, (-np.array(hmax) + np.array(hmax_mc))/np.array(hmax_mc), 'o')
    plt.xlabel("True Energy [TeV]")
    plt.ylabel("(MC hmax - hmax)/(MC hmax)")
    plt.grid()
    plt.semilogx()

'''
    fig=plt.figure()
    ax=plt.gca()
    PlotSummedImages_fromList(ev=9,image=ImageList,mask=MaskList,hillas=HillasList)
    PlotSeveralCams_fromList(ev=9,image=ImageList,hillas=HillasList)

    # check some reconstructed features
    corx = [x['core_x']/u.m for x in ShowerList]
    cory = [x['core_y']/u.m for x in ShowerList]
    alt = [x['alt'][0]/u.deg for x in ShowerList]
    az = [x['az'][0]/u.deg for x in ShowerList]

    # plot azimuth and altitude distribution of the event
    plt.figure()
    plt.hist2d(az,alt,bins=[np.linspace(155,210,50),np.linspace(75,90,50)])
    plt.xlabel("Azimuth angle [deg]")
    plt.ylabel("Altitude angle [deg]")

    # plot core distance
    plt.figure()
    plt.hist2d(corx,cory,bins=[np.linspace(-200,200,50),np.linspace(-200,200,50)])
    plt.show()

'''

