from astropy import units as u
#from astropy.coordinates import Angle
from astropy.time import Time,TimeDelta
from fitsio import FITS      ## change all the code to use this library
from astropy.table import Table
from ctapipe.io.eventsource import EventSource
from ctapipe.io.containers import DataContainer
from ctapipe.instrument import CameraGeometry,OpticsDescription,TelescopeDescription, SubarrayDescription


__all__ = ['HESSfitsIOEventSource']


class HESSfitsIOEventSource(EventSource):
    """
    EventSource for the HESS data in fits file format.

    This class utilises astropy.io.fits to read the HESS fits files, and stores the
    information into the event containers.
    """
    _count = 0

    def __init__(self, config=None, tool=None, **kwargs):
        super().__init__(config=config, tool=tool, **kwargs)
        try:
            import numpy as np
            import os as os
        except ImportError:
            msg = "The `numpy` python module is required to run the pipe"
            self.log.error(msg)
            raise

        self.np = np
        self.os = os

    @staticmethod
    def is_compatible(file_path):
        '''
        Check if the fits file is a valid one and
        if it has the right columns in it
        use build-in commands in astropy.io.fits
        '''
        with FITS(file_path) as hdul:
            try:
                try:
                    hdul['TEVENTS'].has_data() #check if the table has data inside
                except IOError:
                    msg = "Table TEVENTS not found!"
                    #self.log.error(msg)
            except FileNotFoundError:
                msg = "File not found, please provide existing file"
                #self.log.error(msg)
        return True
                

    def __exit__(self, exc_type, exc_val, exc_tb):
        '''
        !!!TO BE DEFINED!!!
        '''
        #hdul.close()
        #del table
        

    def _generator(self):
        with FITS(self.input_url) as hdul:
            # the container is initialized once, and data is replaced within
            # it after each yield
            counter = 0  ## progressive counter
            inc = 0      ## counter that jumps to the the row of the new event
            eventtable = hdul['TEVENTS']
            data = DataContainer()
            data.meta['origin'] = "HESSfits"
            data.meta['input_url'] = self.input_url
            data.meta['max_events'] = self.max_events  #not sure how to handle this if events have a separated entry for each telescope
            obs_id = eventtable.read_header()['OBS_ID']
            reftime = Time(eventtable.read_header()['MJDREFI']+eventtable.read_header()['MJDREFF'],format='mjd')
            while inc < eventtable.get_nrows():
                '''
                make the streaming of events different, it is a bit brutal this way
                '''
                if counter == 0:
                    data.inst.subarray = self._build_subarray_info(obs_id)
                eventstream = eventtable[inc:inc + 4]
                event_id=eventtable['EVENT_ID'][inc]
                data.pointing.azimuth = eventtable['AZ_PNT'][inc][0]  # azimuth per event - needs to be a number
                data.pointing.altitude = eventtable['ALT_PNT'][inc][0]  # altitude per event - needs to be a number
                tels_with_data = eventstream[eventstream['EVENT_ID']==event_id]['TEL_ID']
                data.trig.gps_time = (reftime+TimeDelta(eventtable['TIME'][inc][0], format='sec'))
                data.trig.tels_with_trigger = tels_with_data

                data.count = counter
                data.r0.obs_id = obs_id
                data.r0.event_id = event_id
                data.r0.tels_with_data = tels_with_data
                data.r1.obs_id = obs_id
                data.r1.event_id = event_id
                data.r1.tels_with_data = tels_with_data
                data.dl0.obs_id = obs_id
                data.dl0.event_id = event_id
                data.dl0.tels_with_data = tels_with_data

                # handle telescope filtering by taking the intersection of
                # tels_with_data and allowed_tels
                if len(self.allowed_tels) > 0:
                    selected = tels_with_data & self.allowed_tels
                    if len(selected) == 0:
                        continue  # skip event
                    data.r0.tels_with_data = selected
                    data.r1.tels_with_data = selected
                    data.dl0.tels_with_data = selected

                data.r0.tel.clear()
                data.r1.tel.clear()
                data.dl0.tel.clear()
                data.dl1.tel.clear()

                for tel_id in tels_with_data:
                    ## time of the event for each telescope, copied in the dl0.CameraContainer
                    # data.dl0.tel[tel_id].trigger_time = (TimeDelta(eventstream['TIME'][tel_id - 1], format='sec'))
                    npix = len(data.inst.subarray.tel[tel_id].camera.pix_id)
                    data.dl1.tel[tel_id].image = self.np.zeros(npix)
                    data.dl1.tel[tel_id].image[eventstream['TEL_IMG_IPIX'][tel_id-1]]=eventstream['TEL_IMG_INT'][tel_id-1]
                inc = inc + len(tels_with_data) # update the counter to jump to the next event
                yield data
                counter += 1
        return

    def _build_subarray_info(self, runnumber):
        '''
        !!! TEMPORARY PACTH !!!
        !!! NEEDS IMPROVEMENT !!!
        Needs to open a new fits file because the telescope infos
        are not in the main data file
        TO BE MADE FASTER

        Added also the reading of the camera file "chercam.fits"
        with the position of the pixels

        Parameters
        ----------
        runnumber: run number to match the correct file

        Returns
        -------
        SubarrayDescription :
            instrumental information
        '''
        pathtofile = self.os.path.split(self.input_url)[0]+'/'
        newfilename = pathtofile+"run_00"+str(runnumber)+"_std_zeta_eventlist.fits"  ## FIX THIS!!
        chercamfile = pathtofile+"chercam.fits.gz"
        # pixtab = Table.read(chercamfile,format='ascii') #read the table with the pixel position
        pixtab = Table.read(chercamfile,format='fits') #read the table with the pixel position
        subarray = SubarrayDescription("HESS-I")
        try:
            hdu_array = FITS(newfilename)[2] # open directly the table with the telarray
            teldata = hdu_array.read()
            telescope_ids = list(teldata['TELID'])

            for tel_id in telescope_ids:
                cam=pixtab[(tel_id-1)*960:tel_id*960]
                geom = CameraGeometry(
                                       tel_id,
                                       cam['PIX_ID'],
                                       self.np.array(cam['PIX_POSX'])*u.m,
                                       self.np.array(cam['PIX_POSY'])*u.m,
                                       self.np.array(cam['PIX_AREA'])*u.m*u.m,
                                       pix_type='hexagonal'
                )
                foclen = teldata['FOCLEN'][tel_id-1] * u.m
                mirror_area = 108*u.m*u.m                   # hard coded, NEED FIX! This column is empty in the original fits file!
                num_tiles = 382                             # hard coded, missing in the original file
                optic = OpticsDescription('DC','MST','',foclen,mirror_area,num_tiles)
                
                tel_pos = [
                    teldata['POSX'][tel_id-1],
                    teldata['POSY'][tel_id-1],
                    teldata['POSZ'][tel_id-1]
                ] * u.m
                tel = TelescopeDescription(optic,geom)
                subarray.tels[tel_id] = tel
                subarray.positions[tel_id] = tel_pos
            
            return subarray
        except FileNotFoundError:
            msg = "Secondary file not found, check the presence in same folder!"
            self.log.error(msg)
            raise SystemExit
