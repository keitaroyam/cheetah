/*
Cheetah for EIGER
 Written by Takanori Nakane 

Don't forget to set HDF5_PLUGIN_PATH!
*/

//  Distantly based on cheetah-sacla
//   Created by Anton Barty on 20/1/14.
//   Copyright (c) 2014 Anton Barty. All rights reserved.

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "cheetah.h"

int main(int argc, const char * argv[]) {
	printf("Cheetah for EIGER\n");
	printf(" by Takanori Nakane\n");
	printf(" based on cheetah-sacla by Anton Barty\n");
	printf("\nIf this program fails, make sure HDF5_PLUGIN_PATH points to lz4 plugin!\n");
	
	char filename[1024], cheetahini[1024];
	strcpy(filename, argv[1]);
	strcpy(cheetahini, argv[2]);
   
    printf("Program name: %s\n",argv[0]);
    printf("Input data file: %s\n", filename);
    printf("Cheetah .ini file: %s\n", cheetahini);


	// Parse metadata
	printf("\n** Parsing metadata in the HDF5 file **\n");

	int xpixels = -1, ypixels = -1, beamx = -1, beamy = -1, nimages = -1;
	double pixelsize = -1, wavelength = -1, distance = -1;
	double count_time = -1, frame_time = -1, osc_width = -1;
	char detector_sn[256] = {}, description[256] = {};

	hid_t hdf;	

	hdf = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
	if (hdf < 0) {
		fprintf(stderr, "Failed to open file %s\n", argv[1]);
		return -1;
	}
	printf("Successfully opened %s\n", argv[1]);

	H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/nimages", &nimages);
	printf("Found %d images in the input\n", nimages);

	H5Eset_auto(0, NULL, NULL); // Comment out this line for debugging.

	fprintf(stderr, "Metadata in HDF5:\n");
	H5LTread_dataset_string(hdf, "/entry/instrument/detector/description", description);
	fprintf(stderr, " /entry/instrument/detector/description = %s\n", description);
	H5LTread_dataset_string(hdf, "/entry/instrument/detector/detector_number", detector_sn);
	fprintf(stderr, " /entry/instrument/detector/detector_number = %s\n", detector_sn);
	H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/x_pixels_in_detector", &xpixels);
	H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/y_pixels_in_detector", &ypixels);
	fprintf(stderr, " /entry/instrument/detector/detectorSpecific/{x,y}_pixels_in_detector = (%d, %d) (px)\n",
			xpixels, ypixels);
	H5LTread_dataset_int(hdf, "/entry/instrument/detector/beam_center_x", &beamx);
	H5LTread_dataset_int(hdf, "/entry/instrument/detector/beam_center_y", &beamy);
	fprintf(stderr, " /entry/instrument/detector/beam_center_{x,y} = (%d, %d) (px)\n", beamx, beamy);
	H5LTread_dataset_double(hdf, "/entry/instrument/detector/count_time", &count_time); // in m
	fprintf(stderr, " /entry/instrument/detector/count_time = %f (sec)\n", count_time);
	H5LTread_dataset_double(hdf, "/entry/instrument/detector/frame_time", &frame_time); // in 
	fprintf(stderr, " /entry/instrument/detector/frame_time = %f (sec)\n", frame_time);
	H5LTread_dataset_double(hdf, "/entry/instrument/detector/x_pixel_size", &pixelsize); // in 
	H5LTread_dataset_double(hdf, "/entry/instrument/detector/detector_distance", &distance);
	fprintf(stderr, " /entry/instrument/detector/detector_distance = %f (m)\n", distance);
	fprintf(stderr, " /entry/instrument/detector/x_pixel_size = %f (m)\n", pixelsize);

	// There are many ways to store wavelength!
	H5LTread_dataset_double(hdf, "/entry/instrument/beam/wavelength", &wavelength);
	if (wavelength > 0) {
		fprintf(stderr, " /entry/instrument/beam/wavelength = %f (A)\n", wavelength);
	} else {
		fprintf(stderr, "  /entry/instrument/beam/wavelength not present. Trying another place.\n");
		H5LTread_dataset_double(hdf, "/entry/instrument/monochromator/wavelength", &wavelength);
		if (wavelength > 0) {
			fprintf(stderr, " /entry/instrument/monochromator/wavelength = %f (m)\n", wavelength);
		} else {
			fprintf(stderr, "  /entry/instrument/monochromator/wavelength not present. Trying another place.\n");
			H5LTread_dataset_double(hdf, "/entry/instrument/beam/incident_wavelength", &wavelength);
			if (wavelength > 0) {
				fprintf(stderr, " /entry/instrument/beam/incident_wavelength = %f (m)\n", wavelength);
			}
		}
	}
	if (wavelength < 0) {
		fprintf(stderr, " wavelength not defined!");
	}

	// Same for oscillation range (but not relevant here)
	H5LTread_dataset_double(hdf, "/entry/sample/goniometer/omega_range_average", &osc_width);
	if (osc_width > 0) {
		fprintf(stderr, " /entry/sample/goniometer/omega_range_average = %f (deg)", osc_width);
	} else {
		fprintf(stderr, " oscillation width not defined. \"Start_angle\" is set to 0!\n");
		osc_width = 0;
	}

	
	hid_t entry, group;
	entry = H5Gopen2(hdf, "/entry", H5P_DEFAULT);
	if (entry < 0) {
		fprintf(stderr, "/entry does not exist!\n");
		return -1;
	}
	
	// Check if /entry/data present
	group = H5Gopen2(entry, "data", H5P_DEFAULT);  
	if (group < 0) {
		group = entry; // FIXME: leak!
	}

	int start_block_number = 1;
	if (H5LTfind_dataset(group, "data_000000")) {
		fprintf(stderr, "This dataset starts from data_000000.\n");
		start_block_number = 0;
	} else {
		fprintf(stderr, "This dataset starts from data_000001.\n");
	}

	char data_name[20] = {};
	hid_t data, dataspace;
	int number_per_block = 0;
	
	// Open the first data block to get the number of frames in a block
	snprintf(data_name, 20, "data_%06d", start_block_number); 
	data = H5Dopen2(group, data_name, H5P_DEFAULT);
	dataspace = H5Dget_space(data);
	if (data < 0) {
		fprintf(stderr, "failed to open /entry/%s\n", data_name);
		return -1;
	}
	if (H5Sget_simple_extent_ndims(dataspace) != 3) {
		fprintf(stderr, "Dimension of /entry/%s is not 3!\n", data_name);
		return -1;    
	}
	hsize_t dims[3];
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	number_per_block = dims[0];
	fprintf(stderr, "The number of images per data block is %d.\n", number_per_block);

	H5Sclose(dataspace);
	H5Dclose(data);
	
	//	Initialise Cheetah
	printf("** Setting up Cheetah **\n");
	static uint32_t ntriggers = 0;
	long frameNumber = 0;
    long runNumber = 0;
	static cGlobal cheetahGlobal;
	signed short *buf = (signed short*)malloc(sizeof(signed short) * xpixels * ypixels);

	static time_t startT = 0;
	time(&startT);
    strcpy(cheetahGlobal.configFile, cheetahini);
	cheetahInit(&cheetahGlobal);

    for (frameNumber = 0; frameNumber < nimages; frameNumber++) {
        printf("Processing frame %ld\n", frameNumber);
		
		//  Cheetah: Calculate time beteeen processing of data frames
		time_t	tnow;
		double	dtime, datarate;
		time(&tnow);
		
		dtime = difftime(tnow, cheetahGlobal.tlast);
		if(dtime > 1.) {
			datarate = (frameNumber - cheetahGlobal.lastTimingFrame)/dtime;
			cheetahGlobal.lastTimingFrame = frameNumber;
			time(&cheetahGlobal.tlast);
			cheetahGlobal.datarate = datarate;
		}
			
		// Cheetah: Create a new eventData structure in which to place all information
		cEventData	*eventData;
		eventData = cheetahNewEvent(&cheetahGlobal);
		ntriggers++;
            
        // Cheetah: Populate event structure with meta-data
		eventData->frameNumber = frameNumber;
		eventData->runNumber = runNumber;
		eventData->nPeaks = 0;
		eventData->pumpLaserCode = 0;
		eventData->pumpLaserDelay = 0;
		eventData->photonEnergyeV = 12390 / wavelength;
		eventData->wavelengthA = wavelength;
		eventData->pGlobal = &cheetahGlobal;
		
        // Cheetah: Copy image data into
        
		// 1. open the required data
		int block_number = start_block_number + frameNumber / number_per_block;
		int frame_in_block = frameNumber % number_per_block;
		fprintf(stderr, "frame %d is in block %d, frame %d (all 0-indexed).\n", 
				frameNumber, block_number, frame_in_block);
	
		snprintf(data_name, 20, "data_%06d", block_number); 
		data = H5Dopen2(group, data_name, H5P_DEFAULT);
		dataspace = H5Dget_space(data);
		if (data < 0) {
			fprintf(stderr, "failed to open /entry/%s\n", data_name);
			return -1;
		}
		if (H5Sget_simple_extent_ndims(dataspace) != 3) {
			fprintf(stderr, "Dimension of /entry/%s is not 3!\n", data_name);
			return -1;    
		}
		
		// 2. get the frame
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		hsize_t offset_in[3] = {frame_in_block, 0, 0};
		hsize_t offset_out[3] = {0, 0, 0};
		hsize_t count[3] = {1, ypixels, xpixels};
		hid_t memspace = H5Screate_simple(3, dims, NULL);   
		if (memspace < 0) {
			fprintf(stderr, "failed to create memspace\n");
			return -1;
		}
		int ret;
		ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset_in, NULL, 
								  count, NULL);
		if (ret < 0) {
			fprintf(stderr, "select_hyperslab for file failed\n");
			return -1;
		}
		ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, 
								  count, NULL);
		if (ret < 0) {
			fprintf(stderr, "select_hyperslab for memory failed\n");
			return -1;
		}
		H5Dread(data, H5T_NATIVE_SHORT, memspace, dataspace, H5P_DEFAULT, buf);
		H5Sclose(dataspace);
		H5Sclose(memspace);
		H5Dclose(data);

		// 3. copy to the event data
		long    detID = 0;
		long    pix_nn = cheetahGlobal.detector[detID].pix_nn;
		eventData->detector[detID].data_raw16 = (uint16_t*) calloc(pix_nn, sizeof(uint16_t));

//		int underflow = 0, overflow = 0;
		for(long ii=0; ii<pix_nn; ii++) {
			long tmp = buf[ii];
//			if (ii < 10000) printf("%d ", tmp);
			if (tmp < 0) {
//				underflow++; //printf("%ld ", tmp);
				tmp = 0; 
			} else if (tmp >= SHRT_MAX) { // SHRT_MAX == 32767 seems panel gap. 
//				overflow++;
				tmp = 0; 
			}
			eventData->detector[detID].data_raw16[ii] = (uint16_t)tmp;
		}
//		printf("underflow %d overflow %d\n", underflow, overflow);
		cheetahProcessEventMultithreaded(&cheetahGlobal, eventData);
    }

	// Finish

	cheetahExit(&cheetahGlobal);

	H5Gclose(group);
	H5Fclose(hdf);
	free(buf);
	
	time_t endT;
	time(&endT);
	double dif = difftime(endT,startT);
	std::cout << "time taken: " << dif << " seconds\n";	
	std::cout << "Clean exit!\n";
    return 0;
}

