#
# Combine unit source tsunami initial conditions to make tsunami initial
# conditions for earthquake events.
#
# Here we illustrate stochastic slip, using the S_{NCF} method from:
# Davies et al (2015) Tsunami inundation from heterogeneous earthquake
# slip distributions: Evaluation of synthetic source models. JGR,
# doi:10.1002/2015JB012272
#
# With slight modifications, this code can also generate uniforms slip events
#
#
# Note this code requires that the tutorial in ../source_contours_2_unit_sources
# has been successfully run.
#

library(rptha)

## Input parameters ###########

# Folder containing one directory for each sourcename. Inside the latter
# directories are tif files for each unit source (and no other tif files)
#unit_source_dirname = '../source_contours_2_unit_sources/OUTPUTS/Unit_source_data'
unit_source_dirname = '/home/ignatius.pranantyo/Tsunamis/Learnings/RPTHA_examples/scripts/OUTPUTS/Unit_source_data'

# sourcename. This should be the name of a directory inside
# unit_source_dirname, and also the name of a discretised_source (among
# those contained in all_discretized_source_RDS)
# sourcename = 'alaska'
sourcename = 'jawa__slab2__edited'

# RDS filename containing all discretized source information. The object
# therein should include a list entry corresponding to sourcename
#all_discretized_source_RDS = '../source_contours_2_unit_sources/OUTPUTS/all_discretized_sources.RDS'
all_discretized_source_RDS = '/home/ignatius.pranantyo/Tsunamis/Learnings/RPTHA_examples/scripts/OUTPUTS/all_discretized_sources.RDS'

# Earthquake parameters
desired_Mw = 8.8
target_location = c(110, -9) # Approximate Lon, Lat of rupture (will stochastically vary)
number_of_sffm = 5 # How many stochastic events to make


## end input ##################


discretized_source = readRDS(all_discretized_source_RDS)[[sourcename]]
discretized_source_statistics = discretized_source_approximate_summary_statistics(
    discretized_source)

# Read the raster corresponding to each row in the discretized_source_statistics
unit_source_raster_files = paste(unit_source_dirname, '/', sourcename, '/', 
    sourcename, '_', discretized_source_statistics$downdip_number, '_',
    discretized_source_statistics$alongstrike_number, '.tif', sep="")

if(!(all(file.exists(unit_source_raster_files)))){
    stop('Could not find some unit source raster files')
}

unit_source_rasters = lapply(as.list(unit_source_raster_files), f<-function(x) raster(x))

#
# Make the stochastic slip events. 
# Note many of these options can be varied, e.g. to remove stochasticity in the rupture size or location, or to
# use uniform slip. 
# We often use default values below -- but name them to draw your attention to the existence of these options.
# See ?sffm_make_events_on_discretized_source for further information.
#
stochastic_slip_events = sffm_make_events_on_discretized_source(
    discretized_source_statistics = discretized_source_statistics, 
    target_location = target_location,
    target_event_mw = desired_Mw,
    num_events = number_of_sffm,
    # Enhanced variability in the event location
    vary_peak_slip_location=TRUE, 
    sourcename = sourcename,
    # Use stochastic slip
    uniform_slip = FALSE, 
    # If the desired width is > source-zone width, then in 50% of cases,
    # increase the length to compensate. In the other 50% of cases do not
    # change the length. 
    expand_length_if_width_limited = 'random', 
    # Randomly vary length, width, and corner wavenumber parameters
    use_deterministic_LWkc=FALSE, 
    # If randomly varying length, width, and corner wavenumbers, limit variation to 2SD
    clip_random_parameters_at_2sd = TRUE, 
    # Scaling relation passed to 'Mw_2_rupture_size'
    relation='Strasser', 
    # Do not put the peak slip near the rupture centre (because real
    # earthquakes don't consistently do that)
    peak_slip_location_near_centre=FALSE)

# Store as table
stochastic_slip_events_table = sffm_events_to_table(stochastic_slip_events)


#####################################
### SAVING IMPORTANT OUTPUT FILES ###
#####################################
write.csv(unit_source_raster_files, './OUTPUTS/unit_source_grid_raster_filename_index.csv',
	  row.names=TRUE)
write.csv(stochastic_slip_events_table, './OUTPUTS/stochastic_events_table.csv',
	  row.names=TRUE)



# For each stochastic slip event, compute the vertical deformation by summing
# the unit-sources 
for(i in 1:length(stochastic_slip_events)){

    # Get the contributing event indices and their slip from the table. This data is
    # stored as a seperated string (with separators '-' and '_' respectively),
    # so we need to extract it.
    event_inds = as.numeric(strsplit(stochastic_slip_events_table$event_index_string[i], '-')[[1]])
    event_slip = as.numeric(strsplit(stochastic_slip_events_table$event_slip_string[i], '_')[[1]])

    stopifnot(length(event_inds) == length(event_slip))

    # Sum the raster
    r1 = raster(unit_source_raster_files[1])*0
    for(j in 1:length(event_inds)){
        ei = event_inds[j]
        r1 = r1 + event_slip[j] * raster(unit_source_raster_files[ei])
    }
    # Append to the stochastic_slip_events list
    stochastic_slip_events[[i]]$source_raster = r1
}

# Comment out kajiura, as the repository does not contain elevation data
kajiura = FALSE
if(kajiura){
    # Apply kajiura to the raster. This can be less artefact prone than applying
    # it first to the unit sources, and then summing. Although mathematically the
    # two are equivalent, numerically the unit source approach is somewhat more
    # prone to artefacts

    for(i in 1:length(stochastic_slip_events)){
        stochastic_slip_events[[i]]$source_raster_smooth = 
            kajiura_smooth_raster(
                source_raster=stochastic_slip_events[[i]]$source_raster,
                new_origin=c(209, 58),
                elevation_raster = '../../../../DATA/ELEV/GEBCO_08/gebco_08.nc',
                kj_filter_grid_dxdy = 2000,
                kj_filter_def_threshold=1.0e-02,
                kj_cartesian_buffer = 10000,
                minimum_kj_depth = 10,
                elevation_extraction_x_offset=-360,
                spherical_input = TRUE)
    }
}else{
    # No smoothing, just copy the raster to support the plot below
    for(i in 1:length(stochastic_slip_events)){
        stochastic_slip_events[[i]]$source_raster_smooth = stochastic_slip_events[[i]]$source_raster
    }

}


#
# Make a plot -- this will vary randomly each time the code is run
#
#png('Deformation_plot.png', width=15, height=9, units='in', res=300)
#par(mfrow=c(2,3))
#par(mar=c(3,2,1,1))
#for(i in 1:3) plot(stochastic_slip_events[[i]]$slip_raster, xlab='Along-strike distance (km)', 
#    ylab='Up-dip distance (km)', xlim=c(100, 500))
#for(i in 1:3) {
#    persp(stochastic_slip_events[[i]]$source_raster_smooth, col='skyblue', shade=TRUE, 
#        border=NA, phi=30, theta=-90, scale=FALSE, expand=0.2)
#}
#dev.off()


write_source_rasters_as_tif = TRUE #FALSE
if(write_source_rasters_as_tif){
    #
    # Export each source_raster_smooth as GTiff
    #
    for(i in 1:length(stochastic_slip_events)){
        filename = paste0('./OUTPUTS/source_raster_smooth_', substring(1e+05 + i, 2, 6), '.tif') 
        writeRaster(stochastic_slip_events[[i]]$source_raster_smooth, file=filename, options=c('COMPRESS=DEFLATE'))
    }
}
