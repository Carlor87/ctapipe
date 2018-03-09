from astropy import units as u
from astropy.coordinates import Angle
from astropy.time import Time
from astropy.io import fits
from ctapipe.io.eventsource import EventSource
from ctapipe.io.containers import DataContainer
from ctapipe.instrument import CameraGeometry,OpticsDescription,TelescopeDescription, SubarrayDescription
import gzip
import struct

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
        # NEEDS numpy
        try:
            import numpy as np
        except ImportError:
            msg = "The `numpy` python module is required to run the pipe"
            self.log.error(msg)
            raise

        self.np = np

    @staticmethod
    def is_compatible(file_path):
        '''
        Check if the fits file is a valid one and
        if it has the right columns in it
        use build-in commands in astropy.io.fits
        '''
        with fits.open(file_path,mmap=True) as hdul:
            try:
                try:
                    hdul['TEVENTS'].verify('exception') #should find an exception handling
                except KeyError:
                    msg = "Table TEVENTS not found!"
                    self.log.error(msg)
            except FileNotFoundError:
                msg = "File not found, please provide existing file"
                self.log.error(msg)
        return True;
                

    def __exit__(self, exc_type, exc_val, exc_tb):
        '''
        !!!TO BE BETTER DEFINED!!!
        needs both to close the file
        and to send the table to garbage collection.
        Needed because file opened with mmap=True
        See astropy.io.fits documentation
        '''
        #hdul.close()
        #del table
        

    def _generator(self):
        with fits.open(self.input_url,mmap=True,mode='denywrite') as hdul:
            # the container is initialized once, and data is replaced within
            # it after each yield
            counter = 0  ## progressive counter
            inc = 0      ## counter that jumps to the the row of the new event
            eventtable = hdul['TEVENTS']
            eventstream = eventtable.data
            data = DataContainer()
            data.meta['origin'] = "HESSfits"

            data.meta['input_url'] = self.input_url
            data.meta['max_events'] = self.max_events  #not sure how to handle this if events have a separated entry for each telescope
            obs_id = eventtable.header['OBS_ID']
            
            data.inst.subarray = self._build_subarray_info(obs_id)
            data.pointing.azimuth = eventtable.header['AZ_PNT']    #average azimuth for run 
            data.pointing.altitude = eventtable.header['ALT_PNT']  #average zenith for run (ideally per event)

            while inc < (len(eventstream['EVENT_ID'])):
                '''
                make the streaming of events different, it is a bit brutal this way
                '''
                event_id=eventstream['EVENT_ID'][inc] 
                tels_with_data = eventstream[eventstream['EVENT_ID']==event_id].TEL_ID
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

                telc=0 #counter for the number of telescopes
                for tel_id in tels_with_data:
                    npix = len(data.inst.subarray.tel[tel_id].camera.pix_id)
                    countcam = inc+telc #ausiliar counter for the correct row of the image
                    data.dl1.tel[tel_id].image = self.np.zeros(npix)
                    data.dl1.tel[tel_id].image[eventstream['TEL_IMG_IPIX'][countcam]]=eventstream['TEL_IMG_INT'][countcam]
                    telc+=1
                
                inc = inc + len(tels_with_data) # update the counter to jump to the next event
                yield data
                counter += 1
                

        return

    def _build_subarray_info(self, runnumber):
        """
        !!! TEMPORARY PACTH !!!
        !!! NEEDS IMPROVEMENT !!!
        Needs to open a new fits file because the telescope infos
        are not in the main data file

        Parameters
        ----------
        runnumber: run number to match the correct file

        Returns
        -------
        SubarrayDescription :
            instrumental information
        """
        newfilename = "run_00"+str(runnumber)+"_std_zeta_eventlist.fits"  ##FIX THIS!!
        subarray = SubarrayDescription("HESS-I")
        try:
            hdu_array = fits.open(newfilename)[2] #open directly the table with the telarray    
            teldata = hdu_array.data
            telescope_ids = list(teldata['TELID'])

            for tel_id in telescope_ids:
                geom = CameraGeometry.from_name("HESS-I")   # import HESS-I camera geometry
                foclen = teldata['FOCLEN'][tel_id-1] * u.m
                mirror_area = 108*u.m*u.m                   # hard coded, NEED FIX! This column is empty in the original fits file!
                num_tiles = 382                             # hard coded, missing in the original file
                optic = OpticsDescription('DC','MST','',foclen,mirror_area,num_tiles)
                
                tel_pos = [teldata['POSX'][tel_id-1] * u.m,teldata['POSY'][tel_id-1] * u.m]
                tel = TelescopeDescription(optic,geom)
                subarray.tels[tel_id] = tel
                subarray.positions[tel_id] = tel_pos
            
            return subarray
        except FileNotFoundError:
            msg = "Secondary file not found, check the presence in same folder!"
            self.log.error(msg)
            raise SystemExit
