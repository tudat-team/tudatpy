
#include <string>


namespace tudatpy {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else {
        return "No documentation found.";
    }

}
//
//
//
//namespace interface {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//
//namespace spice {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//    } else if(name == "SpiceEphemeris") {
//         return R"(
//
//        Ephemeris derived class which retrieves the state of a body directly from the SPICE library.
//
//        Ephemeris derived class which retrieves the state of a body directly from the SPICE library.
//        The body of which the ephemeris is to be retrieved, as well as the origin and orientation
//        of the reference frame in which the states are returned, and any corrections that are
//        applied, are defined once during object construction.
//
//     )";
//
//
//    } else if(name == "SpiceEphemeris.__init__" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Constructor, sets the input variables for the calls to the spice function to retrieve state.
//
//        Parameters
//        ----------
//        target_body_name
//            Name of body of which the ephemeris is to be calculated.
//        observer_body_name
//            Name of body relative to which the ephemeris is to be calculated.
//        correct_for_stellar_aberration
//            Boolean whether to correct for stellar Aberration in retrieved values of (observed state).
//
//        correct_for_light_time_aberration
//            Boolean whether to correct for light time in retrieved values of (observed state).
//
//        converge_ligh_time_aberration
//            Boolean whether to use single iteration or max. 3 iterations for calculating light time.
//
//        reference_frame_name
//            Name of the reference frame in which the epehemeris is to be calculated.
//
//        reference_julian_day
//            Reference julian day w.r.t. which ephemeris is evaluated.
//
//    )";
//
//
//
//    } else if(name == "SpiceEphemeris.get_cartesian_state" && variant==0) {
//            return R"(
//
//        Get Cartesian state from ephemeris.
//
//         Returns Cartesian state from ephemeris at given Julian day.
//
//        Parameters
//        ----------
//        seconds_since_epoch : float
//            Seconds since epoch at which ephemeris is to be evaluated.
//    )";
//
//
//
//
//
//    } else if(name == "convert_julian_date_to_ephemeris_time" && variant==0) {
//            return R"(
//
//        Convert a Julian date to ephemeris time (equivalent to TDB in Spice).
//
//        The following math is for documentation demonstration purposes
//
//        .. math:: X(e^{j\omega } ) = x(n)e^{ - j\omega n}
//
//        \f$ f(x) = a + b \f$
//
//        Function to convert a Julian date to ephemeris time, which is
//        equivalent to barycentric dynamical time. A leap second kernel
//        must have been loaded to use this function.
//
//
//        Parameters
//        ----------
//        julian_date : int
//            Julian date that is to be converted to ephemeris time.
//
//        Returns
//        -------
//        ephemeris_time : float    Julian date calculated from ephemeris time.
//
//    )";
//
//
//
//    } else if(name == "convert_ephemeris_time_to_julian_date" && variant==0) {
//            return R"(
//
//        Convert ephemeris time (equivalent to TDB) to a Julian date.
//
//        Function to convert ephemeris time, which is nearly equal to
//        barycentric dynamical time, to the Julian date. A leap second
//        kernel must have been loaded to use this function.
//
//
//        Parameters
//        ----------
//        ephemeris_time : float
//            Ephemeris time that is to be converted to Julian date.
//
//        Returns
//        -------
//        julian_date : float    Julian date calculated from ephemeris time.
//
//    )";
//
//
//
//    } else if(name == "convert_date_string_to_ephemeris_time" && variant==0) {
//            return R"(
//
//        Converts a date string to ephemeris time.
//
//        Function to convert a date string, for instance
//        1988 June 13, 3:29:48 to ephemeris time, wrapper for `str2et_c`
//        spice function.
//
//
//        Parameters
//        ----------
//        date_string : str
//            String representing the date. See documentation of spice
//            function `str2et_c` for details on supported formats.
//
//
//        Returns
//        -------
//        ephemeris_time : str    Ephemeris time corresponding to given date_string.
//
//    )";
//
//
//
//    } else if(name == "get_body_cartesian_state_at_epoch" && variant==0) {
//            return R"(
//
//        Get Cartesian state of a body, as observed from another body.
//
//        This function returns the state of a body, relative to another
//        body, in a frame specified by the user. Corrections for light-time
//        correction and stellar aberration can be applied to obtain the
//        state of one of the bodies, as observed from the other. Wrapper
//        for `spkezr_c` spice function.
//
//
//        Parameters
//        ----------
//        target_body_name : str
//            Name of the body of which the state is to be obtained. A kernel
//            with the ephemeris of this body must have been loaded. The
//            string must be a spice-recognized name or ID.
//
//        observer_body_name : str
//            Name of the body relative to which the state is to be obtained.
//            A kernel with the ephemeris of this body must have been loaded.
//            The string must be a spice-recognized name or ID.
//
//        reference_frame_name : str
//            The spice-recognized name of the reference frame in which the
//            state is to be returned. Spice kernel(s) required to perform
//            the necessary conversion from the states of the target and
//            observer bodies to this frame need to have been loaded.
//
//        aberration_corrections : str
//            Setting for correction for setting corrections. See Spice
//            documentation for extended discussion.
//            Short summary:
//
//            - NONE: none
//            - LT: light time corrected (one iteration for calculation)
//            - CN: light time corrected (multiple iterations, max 3) for calculation
//            - S: Stellar aberration corrected.
//            - XLT and XCN: can be provided to make the ephemeris time input argument the transmission time, instead of reception time. Arguments can be combined (i.e."LT+S" or "XCN+S").
//
//        ephemeris_time : float
//            Observation time (or transmission time of observed light, see description
//            of aberrationCorrections).
//
//
//        Returns
//        -------
//        cartesian_state_vector : np.ndarray[6,]    Cartesian state vector (x,y,z, position+velocity).
//
//    )";
//
//
//
//    } else if(name == "get_body_cartesian_position_at_epoch" && variant==0) {
//            return R"(
//
//        Get Cartesian position of a body, as observed from another body.
//
//        This function returns the position of a body, relative to another
//        body, in a frame specified by the user. Corrections for light-time
//        correction and stellar aberration can be applied to obtain the
//        state of one of the bodies, as observed from the other. Wrapper
//        for `spkpos_c` spice function.
//
//
//        Parameters
//        ----------
//        target_body_name : str
//            Name of the body of which the state is to be obtained. A kernel
//            with the ephemeris of this body must have been loaded. The
//            string must be a spice-recognized name or ID.
//
//        observer_body_name : str
//            Name of the body relative to which the state is to be obtained.
//            A kernel with the ephemeris of this body must have been loaded.
//            The string must be a spice-recognized name or ID.
//
//        reference_frame_name : str
//            The spice-recognized name of the reference frame in which the
//            state is to be returned. Spice kernel(s) required to perform
//            the necessary conversion from the states of the target and
//            observer bodies to this frame need to have been loaded.
//
//        aberration_corrections : str
//            Setting for correction for setting corrections. See Spice
//            documentation for extended discussion.
//            Short summary:
//
//            - NONE: none
//            - LT: light time corrected (one iteration for calculation)
//            - CN: light time corrected (multiple iterations, max 3) for calculation,
//            - S: Stellar aberration corrected.
//            - XLT and XCN: can be provided to make the ephemeris time input argument the transmission time, instead of reception time. Arguments can be combined (i.e."LT+S" or "XCN+S").
//
//        ephemeris_time : float
//            Observation time (or transmission time of observed light, see description
//            of aberrationCorrections).
//
//    )";
//
//
//
//    } else if(name == "get_cartesian_state_from_tle_at_epoch" && variant==0) {
//            return R"(
//
//        Get Cartesian state of a satellite from its two-line element set at a specified epoch.
//
//        This function retrieves the state of a satellite at a certain epoch
//        by propagating the SGP or SDP models (near-Earth resp. deep space)
//        with the given two-line elements (TLE). This function serves as a
//        wrapper for the `ev2lin_` function in CSpice.
//
//
//        Parameters
//        ----------
//        epoch : float
//            Time in seconds since J2000 at which the state is to be retrieved.
//        tle : :class:`~tudatpy.kernel.astro.ephemerides.Tle`
//            Shared pointer to a Tle object containing the SGP/SDP model parameters as derived from the element set.
//
//        Returns
//        -------
//        cartesian_state_vector : np.ndarray[6,]    Cartesian state vector (x,y,z, position+velocity).
//
//    )";
//
//
//
//    } else if(name == "compute_rotation_quaternion_between_frames" && variant==0) {
//            return R"(
//
//        Compute quaternion of rotation between two frames.
//
//        This function computes the quaternion of rotation between two
//        frames at a given time instant. kernels defining the two frames,
//        as well as any required intermediate frames, at the requested
//        time must have been loaded. Wrapper for `pxform_c` spice function.
//
//
//        Parameters
//        ----------
//        original_frame
//            Reference frame from which the rotation is made.
//        new_frame
//            Reference frame to which the rotation is made.
//        ephemeris_time
//            Value of ephemeris time at which rotation is to be determined.
//
//        Returns
//        -------
//        Rotation quaternion from original to new frame at given time.
//
//    )";
//
//
//
//    } else if(name == "compute_rotation_matrix_derivative_between_frames" && variant==0) {
//            return R"(
//
//        Computes time derivative of rotation matrix between two frames.
//
//        This function computes the derivative of the rotation matrix
//        between two frames at a given time instant. kernels defining the
//        two frames, as well as any required intermediate frames, at the
//        requested time must have been loaded. Wrapper for (part of) `sxform_c` spice function.
//
//
//        Parameters
//        ----------
//        original_frame
//            Reference frame from which the rotation is made.
//        new_frame
//            Reference frame to which the rotation is made.
//        ephemeris_time
//            Value of ephemeris time at which rotation is to be determined.
//
//        Returns
//        -------
//        Time derivative of rotation matrix from original to new frame at given time.
//
//    )";
//
//
//
//    } else if(name == "get_angular_velocity_vector_of_frame_in_original_frame" && variant==0) {
//            return R"(
//
//        Computes the angular velocity of one frame w.r.t. to another frame.
//
//        Computes the angular velocity of one frame w.r.t. to another frame.
//        at a given time instant. kernels defining the two frames, as well
//        as any required intermediate frames, at the requested time must
//        have been loaded. Wrapper for `xf2rav_c`_ spice function (utilizing `sxform_c`_).
//
//        .. _`xf2rav_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/xf2rav_c.html
//        .. _`sxform_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sxform_c.html
//
//
//        Parameters
//        ----------
//        original_frame
//            Reference frame from which the rotation is made.
//        new_frame
//            Reference frame to which the rotation is made.
//        ephemeris_time
//            Value of ephemeris time at which rotation is to be determined.
//
//        Returns
//        -------
//        Angular velocity of newFrame w.r.t. originalFrame, expressed in originalFrame.
//
//    )";
//
//
//
//    } else if(name == "get_body_properties" && variant==0) {
//            return R"(
//
//        Get property of a body from Spice.
//
//        Function to retrieve a property of a body from Spice, wraps the bodvrd_c Spice function.
//
//
//        Parameters
//        ----------
//        body_name
//            Name of the body of which the property is to be retrieved.
//        property
//            Name of the property that is to be retrieved. Naming conventions can be found
//            in the `bodvrd_c`_ function documentation.
//
//            .. _`bodvrd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html
//
//        maximum_number_of_values : int
//            Number of values by which the property is expressed (i.e. 1 for
//            gravitational parameter, 3 for tri-axial ellipsoid principal axes).
//
//
//        Returns
//        -------
//        Property value(s) expressed in an STL vector of doubles.
//
//    )";
//
//
//
//    } else if(name == "get_body_gravitational_parameter" && variant==0) {
//            return R"(
//
//        Get gravitational parameter of a body.
//
//        This function retrieves the gravitational parameter of a body.
//        Wraps the `bodvrd_c`_ spice function with "GM" as property type.
//
//        .. _`bodvrd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html
//
//
//        Parameters
//        ----------
//        body
//            Name of the body of which the parameter is to be retrieved.
//
//        Returns
//        -------
//        Gravitational parameter of requested body.
//
//    )";
//
//
//
//    } else if(name == "get_average_radius" && variant==0) {
//            return R"(
//
//        Get the (arithmetic) mean of the three principal axes of the tri-axial ellipsoid shape.
//
//        Returns the (arithmetic) mean of the three principal axes of the
//        tri-axial ellipsoid shape of the requested body. Uses the `bodvrd_c` spice function with "RADII" as property type.
//
//        .. _`bodvrd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html
//
//
//        Parameters
//        ----------
//        body
//            Name of the body of which the average radius is to be retrieved.
//
//        Returns
//        -------
//        Arithmetic mean of principal axes of tri-axial ellipsoid shape model of body.
//
//    )";
//
//
//
//    } else if(name == "convert_body_name_to_naif_id" && variant==0) {
//            return R"(
//
//        Convert a body name to its NAIF identification number.
//
//        This function converts a body name to its NAIF identification
//        number. The NAIF id number is required for a number of spice
//        functions, whereas the name is easily interpretable by the user.
//        Wrapper for the ``bods2c_c`` function.
//
//
//        Parameters
//        ----------
//        body_name
//            Name of the body for which NAIF id is to be retrieved.
//
//        Returns
//        -------
//        NAIF id number for the body with bodyName.
//
//    )";
//
//
//
//    } else if(name == "check_body_property_in_kernel_pool" && variant==0) {
//            return R"(
//
//        Check if a certain property of a body is in the kernel pool.
//
//        This function checks if a certain property of a body is in the
//        kernel pool. These properties are defined in PCK kernels. Their
//        names are given in the kernel file, typical names can be found in
//        the Spice documentation. Wrapper for the `bodfnd_c`_ function.
//
//        .. _`bodfnd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodfnd_c.html
//
//
//        Parameters
//        ----------
//        body_name
//            Name of the body of which the property is to be checked.
//        body_property
//            Name of the property of which the presence is to be checked, not case-sensitive.
//
//        Returns
//        -------
//        bool
//            True if property is in pool, false if not.
//
//    )";
//
//
//
//    } else if(name == "get_standard_kernels" && variant==0) {
//            return R"(
//
//        Get the paths to the default legacy kernels.
//
//    )";
//
//
//
//    } else if(name == "load_standard_kernels" && variant==0) {
//            return R"(
//
//        Load the default legacy kernels.
//
//
//        Parameters
//        ----------
//        kernel_paths : List[str]
//            Optional addition kernels to be loaded.
//    )";
//
//
//
//    } else if(name == "get_total_count_of_kernels_loaded" && variant==0) {
//            return R"(
//
//        Get the number of spice kernels currently loaded.
//
//        This function returns the amount of Spice kernels that are loaded
//        into the kernel pool. The same kernel can be loaded multiple times.
//        Wrapper for the `ktotal_c`_ function.
//
//        .. _`ktotal_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ktotal_c.html
//
//
//        Returns
//        -------
//        n_kernels : int    Number of spice kernels currently loaded.
//
//    )";
//
//
//
//    } else if(name == "load_kernel" && variant==0) {
//            return R"(
//
//        Loads a Spice kernel into the pool.
//
//        This function loads a Spice kernel into the kernel pool, from which
//        it can be used by the various internal spice routines. Matters
//        regarding the manner in which Spice handles different kernels
//        containing the same information can be found in the spice required
//        reading documentation, kernel section. Wrapper for the `furnsh_c`_
//        function.
//
//        .. _`furnsh_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/furnsh_c.html
//
//
//        Parameters
//        ----------
//        file_path : str
//            Path to the spice kernel to be loaded.
//    )";
//
//
//
//    } else if(name == "clear_kernels" && variant==0) {
//            return R"(
//
//        Clear all loaded spice kernels.
//
//        This function removes all Spice kernels from the kernel pool.
//        Wrapper for the `kclear_c`_ function.
//
//        .. _`kclear_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/kclear_c.html
//
//
//        Returns
//        -------
//        None
//            None
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//
//
//
//}
//
//
//
//
//
//namespace simulation {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//
//namespace environment_setup {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//
//namespace ephemeris {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//    } else if(name == "BodiesWithEphemerisData") {
//         return R"(
//
//        Enumeration of bodies with ephemeris data.
//
//        Enumeration of bodies with ephemeris data.
//
//     )";
//
//
//
//    } else if(name == "EphemerisSettings") {
//         return R"(
//
//        Base class for providing settings for ephemeris model.
//
//        Functional (base) class for settings of ephemeris models that require no information in addition to their type (and frame origin and orientation).
//        Ephemeris model classes requiring additional information must be created using an object derived from this class.
//
//     )";
//
//
//
//    } else if(name == "ScaledEphemerisSettings") {
//         return R"(
//
//        Class for defining settings from scaling existing ephemeris settings.
//
//        `EphemerisSettings` derived class for a new ephemeris created from scaling an existing ephemeris settings object. It allows the user to apply a scaling factor to the resulting Cartesian states (for instance for an uncertainty analysis).
//     )";
//
//
//
//    } else if(name == "DirectSpiceEphemerisSettings") {
//         return R"(
//
//        Class for defining settings of an ephemeris linked directly to Spice.
//
//        `EphemerisSettings` derived class for ephemeris which are directly linked to Spice.
//     )";
//
//
//
//    } else if(name == "InterpolatedSpiceEphemerisSettings") {
//         return R"(
//
//        Class for defining settings of an ephemeris interpolated from Spice data.
//
//        `DirectSpiceEphemerisSettings` derived class for setting ephemerides to be created from interpolated Spice ephemeris data.
//     )";
//
//
//
//    } else if(name == "ApproximatePlanetPositionSettings") {
//         return R"(
//
//        Class for creating settings of approximate ephemeris for major planets.
//
//        `EphemerisSettings` derived class for approximate ephemeris for major planets as inplemented in ApproximatePlanetPositions class and derived class (described on http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf).
//     )";
//
//
//
//    } else if(name == "ConstantEphemerisSettings") {
//         return R"(
//
//        Class for defining settings of constant ephemerides.
//
//        `EphemerisSettings` derived class for ephemerides producing a constant (time-independent) state.
//     )";
//
//
//
//    } else if(name == "CustomEphemerisSettings") {
//         return R"(
//
//        Class for defining settings of a custom ephemeris.
//
//        `EphemerisSettings` derived class for ephemerides which represent an ideal Kepler orbit.
//     )";
//
//
//
//    } else if(name == "TabulatedEphemerisSettings") {
//         return R"(
//
//        Class for defining settings of ephemeris to be created from tabulated data.
//
//        `EphemerisSettings` derived class for ephemeris created from tabulated data. The provided data is interpolated into ephemerides.
//     )";
//
//
//
//
//    } else if(name == "direct_spice" && variant==0) {
//            return R"(
//
//        Factory function for creating ephemeris model settings entirely from Spice.
//
//        Factory function for settings object, defining ephemeris model directly and entirely from Spice.
//        Requires an appropriate Spice kernel to be loaded.
//        This function creates an instance of an `EphemerisSettings` derived `DirectSpiceEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        frame_origin : str, default='SSB'
//            Origin of frame in which ephemeris data is defined.
//        frame_orientation : str, default='ECLIPJ2000'
//            Orientation of frame in which ephemeris data is defined.
//        body_name_to_use : str, default = ""
//            ?
//
//        Returns
//        -------
//        DirectSpiceEphemerisSettings
//            None
//
//    )";
//
//
//
//    } else if(name == "interpolated_spice" && variant==0) {
//            return R"(
//
//        Factory function for creating ephemeris model settings using interpolated Spice data.
//
//
//        Parameters
//        ----------
//        initial_time : float
//            Initial time from which interpolated data from Spice should be created.
//        final_time : float
//            Final time from which interpolated data from Spice should be created.
//        time_step : float
//            Time step with which interpolated data from Spice should be created.
//        frame_origin : str, default='SSB'
//            Origin of frame in which ephemeris data is defined.
//        frame_orientation : str, default='ECLIPJ2000'
//            Orientation of frame in which ephemeris data is defined.
//        interpolator_settings : std::make_shared< interpolators::InterpolatorSettings >, default=std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 )
//            Settings to be used for the state interpolation.
//        body_name_to_use : str, default = ""
//            ?
//
//        Returns
//        -------
//        InterpolatedSpiceEphemerisSettings
//            None
//
//    )";
//
//
//
//    } else if(name == "approximate_planet_positions" && variant==0) {
//            return R"(
//
//        Factory function for creating approximate ephemeris model settings for major planets.
//
//        Factory function for settings object, defining approximate ephemeris model for major planets.
//        In this highly simplified ephemeris model, Keplerian elements of the major solar system bodies are modelled as linear functions of time and several sinusoidal variations (described on http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf).
//        #(%! only given in ctor, not FF) Three-dimensional and circular-coplanar approximation may be used.
//        This function creates an instance of an `EphemerisSettings` derived `ApproximatePlanetPositionsSettings` object.
//
//
//        Parameters
//        ----------
//        body_name_to_use : str
//            String that is attempted to be matched to an identifier for the body that the ephemeris is to be created for.
//
//        Returns
//        -------
//        ApproximatePlanetPositionSettings
//            None
//
//    )";
//
//
//    } else if(name == "approximate_planet_positions" && variant==1) {
//            return R"(
//
//        Factory function for creating approximate ephemeris model settings for major planets.
//
//        Factory function for settings object, defining approximate ephemeris model for major planets.
//        In this highly simplified ephemeris model, Keplerian elements of the major solar system bodies are modelled as linear functions of time and several sinusoidal variations (described on http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf).
//        #(%! only given in ctor, not FF) Three-dimensional and circular-coplanar approximation may be used.
//        This function creates an instance of an `EphemerisSettings` derived `ApproximatePlanetPositionsSettings` object.
//
//
//        Parameters
//        ----------
//        None
//            None
//
//        Returns
//        -------
//        ApproximatePlanetPositionSettings
//            None
//
//    )";
//
//
//
//    } else if(name == "approximate_planet_positions" && variant==0) {
//            return R"(
//
//        Factory function for creating approximate ephemeris model settings for major planets.
//
//        Factory function for settings object, defining approximate ephemeris model for major planets.
//        In this highly simplified ephemeris model, Keplerian elements of the major solar system bodies are modelled as linear functions of time and several sinusoidal variations (described on http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf).
//        #(%! only given in ctor, not FF) Three-dimensional and circular-coplanar approximation may be used.
//        This function creates an instance of an `EphemerisSettings` derived `ApproximatePlanetPositionsSettings` object.
//
//
//        Parameters
//        ----------
//        body_name_to_use : str
//            String that is attempted to be matched to an identifier for the body that the ephemeris is to be created for.
//
//        Returns
//        -------
//        ApproximatePlanetPositionSettings
//            None
//
//    )";
//
//
//    } else if(name == "approximate_planet_positions" && variant==1) {
//            return R"(
//
//        Factory function for creating approximate ephemeris model settings for major planets.
//
//        Factory function for settings object, defining approximate ephemeris model for major planets.
//        In this highly simplified ephemeris model, Keplerian elements of the major solar system bodies are modelled as linear functions of time and several sinusoidal variations (described on http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf).
//        #(%! only given in ctor, not FF) Three-dimensional and circular-coplanar approximation may be used.
//        This function creates an instance of an `EphemerisSettings` derived `ApproximatePlanetPositionsSettings` object.
//
//
//        Parameters
//        ----------
//        None
//            None
//
//        Returns
//        -------
//        ApproximatePlanetPositionSettings
//            None
//
//    )";
//
//
//
//    } else if(name == "constant" && variant==0) {
//            return R"(
//
//        Factory function for creating constant ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model with a constant, time-independent state.
//        This function creates an instance of an `EphemerisSettings` derived `constantEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        constant_state : numpy.ndarray
//            Constant state that will be provided as output of the ephemeris at all times.
//        frame_origin : str, default='SSB'
//            Origin of frame in which ephemeris data is defined.
//        frame_orientation : str, default='ECLIPJ2000'
//            Orientation of frame in which ephemeris data is defined.
//
//        Returns
//        -------
//        ConstantEphemerisSettings
//
//
//    )";
//
//
//
//    } else if(name == "custom" && variant==0) {
//            return R"(
//
//        Factory function for creating custom ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model with a custom state.
//        This allows the user to provide an custom state function as ephemeris model.
//        The state function (pointer) must be taking a time (float) as input and returning the Cartesian state (numpy.ndarray).
//        This function creates an instance of an `EphemerisSettings` derived `customEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        custom_state_function
//            Function returning the state as a function of time.
//        frame_origin : str, default='SSB'
//            Origin of frame in which ephemeris data is defined.
//        frame_orientation : str, default='ECLIPJ2000'
//            Orientation of frame in which ephemeris data is defined.
//
//        Returns
//        -------
//        CustomEphemerisSettings
//
//
//    )";
//
//
//
//    } else if(name == "keplerian" && variant==0) {
//            return R"(
//
//        Factory function for creating Keplerian ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model which represents an ideal Kepler orbit from the given Kepler elements.
//        These are taken as the elements at the ``initial_state_epoch`` and propagated to any other time using the provided ``central_body_gravitational_parameter``.
//        See Frame/State Transformations (`link`) for more details on orbital elements in Tudat.
//        This function creates an instance of an `EphemerisSettings` derived `KeplerEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        initial_state_in_keplerian_elements : numpy.ndarray
//            Kepler elements at time epochOfInitialState.
//        initial_state_epoch : float
//            Time at which initialStateInKeplerianElements represents the Keplerian state.
//        central_body_gravitational_parameter : float
//            Gravitational parameter of the central body that is used in the computations.
//        frame_origin : str, default='SSB'
//            Origin of frame in which ephemeris data is defined.
//        frame_orientation : str, default='ECLIPJ2000'
//            Orientation of frame in which ephemeris data is defined.
//        root_finder_absolute_tolerance : float
//            Convergence tolerance for root finder used to convert mean to eccentric anomaly on each call to getCartesianState.
//        root_finder_maximum_number_of_iterations : float
//            Maximum iteration for root finder used to convert mean to eccentric anomaly on each call to getCartesianState.
//
//        Returns
//        -------
//        KeplerEphemerisSettings
//
//
//    )";
//
//
//
//    } else if(name == "keplerian_from_spice" && variant==0) {
//            return R"(
//
//        Factory function for creating Keplerian ephemeris model settings with initial state from Spice.
//
//        Factory function for settings object, defining ephemeris model which represents an ideal Kepler orbit from an initial state from Spice.
//        The Kepler elements inferred from the initial state are propagated to any other time using the provided ``central_body_gravitational_parameter``.
//        See Frame/State Transformations (`link`) for more details on orbital elements in Tudat.
//        This function creates an instance of an `EphemerisSettings` derived `KeplerEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        body : str
//            Name of body for which to create ephemeris settings and infer initial state from Spice.
//        initial_state_epoch : float
//            Time at which initialStateInKeplerianElements represents the Keplerian state.
//        central_body_gravitational_parameter : float
//            Gravitational parameter of the central body that is used in the computations.
//        frame_origin : str, default='SSB'
//            Origin of frame in which ephemeris data is defined.
//        frame_orientation : str, default='ECLIPJ2000'
//            Orientation of frame in which ephemeris data is defined.
//        root_finder_absolute_tolerance : float
//            Convergence tolerance for root finder used to convert mean to eccentric anomaly on each call to getCartesianState.
//        root_finder_maximum_number_of_iterations : float
//            Maximum iteration for root finder used to convert mean to eccentric anomaly on each call to getCartesianState.
//
//        Returns
//        -------
//        KeplerEphemerisSettings
//
//
//    )";
//
//
//
//    } else if(name == "scaled" && variant==0) {
//            return R"(
//
//        Factory function for creating scaled ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
//        The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
//        This function creates an instance of an `EphemerisSettings` derived `ScaledEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        unscaled_ephemeris_settings : EphemerisSettings
//            Sets base settings of ephemeris to be scaled.
//        scaling_constant : float
//            Constant scaling factor to be applied to all elements of the Cartesian state.
//        is_scaling_absolute : bool, default=false
//            Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
//
//        Returns
//        -------
//        ScaledEphemerisSettings
//
//
//    )";
//
//
//    } else if(name == "scaled" && variant==1) {
//            return R"(
//
//        Factory function for creating scaled ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
//        The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
//        This function creates an instance of an `EphemerisSettings` derived `ScaledEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        unscaled_ephemeris_settings : EphemerisSettings
//            Sets base settings of ephemeris to be scaled.
//        scaling_vector : numpy.ndarray
//            Vector containing scaling factors to be applied to each element of the Cartesian state.
//        is_scaling_absolute : bool, default=false
//            Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
//
//        Returns
//        -------
//        ScaledEphemerisSettings
//
//
//    )";
//
//
//    } else if(name == "scaled" && variant==2) {
//            return R"(
//
//        Factory function for creating scaled ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
//        The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
//        This function creates an instance of an `EphemerisSettings` derived `ScaledEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        unscaled_ephemeris_settings : EphemerisSettings
//            Sets base settings of ephemeris to be scaled.
//        scaling_vector_function : function < numpy.ndarray >
//            Function returning a vector with the scaling factors to be applied to each element of the Cartesian state.
//        is_scaling_absolute : bool, default=false
//            Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
//
//        Returns
//        -------
//        ScaledEphemerisSettings
//
//
//    )";
//
//
//
//    } else if(name == "scaled" && variant==0) {
//            return R"(
//
//        Factory function for creating scaled ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
//        The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
//        This function creates an instance of an `EphemerisSettings` derived `ScaledEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        unscaled_ephemeris_settings : EphemerisSettings
//            Sets base settings of ephemeris to be scaled.
//        scaling_constant : float
//            Constant scaling factor to be applied to all elements of the Cartesian state.
//        is_scaling_absolute : bool, default=false
//            Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
//
//        Returns
//        -------
//        ScaledEphemerisSettings
//
//
//    )";
//
//
//    } else if(name == "scaled" && variant==1) {
//            return R"(
//
//        Factory function for creating scaled ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
//        The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
//        This function creates an instance of an `EphemerisSettings` derived `ScaledEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        unscaled_ephemeris_settings : EphemerisSettings
//            Sets base settings of ephemeris to be scaled.
//        scaling_vector : numpy.ndarray
//            Vector containing scaling factors to be applied to each element of the Cartesian state.
//        is_scaling_absolute : bool, default=false
//            Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
//
//        Returns
//        -------
//        ScaledEphemerisSettings
//
//
//    )";
//
//
//    } else if(name == "scaled" && variant==2) {
//            return R"(
//
//        Factory function for creating scaled ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
//        The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
//        This function creates an instance of an `EphemerisSettings` derived `ScaledEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        unscaled_ephemeris_settings : EphemerisSettings
//            Sets base settings of ephemeris to be scaled.
//        scaling_vector_function : function < numpy.ndarray >
//            Function returning a vector with the scaling factors to be applied to each element of the Cartesian state.
//        is_scaling_absolute : bool, default=false
//            Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
//
//        Returns
//        -------
//        ScaledEphemerisSettings
//
//
//    )";
//
//
//
//    } else if(name == "scaled" && variant==0) {
//            return R"(
//
//        Factory function for creating scaled ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
//        The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
//        This function creates an instance of an `EphemerisSettings` derived `ScaledEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        unscaled_ephemeris_settings : EphemerisSettings
//            Sets base settings of ephemeris to be scaled.
//        scaling_constant : float
//            Constant scaling factor to be applied to all elements of the Cartesian state.
//        is_scaling_absolute : bool, default=false
//            Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
//
//        Returns
//        -------
//        ScaledEphemerisSettings
//
//
//    )";
//
//
//    } else if(name == "scaled" && variant==1) {
//            return R"(
//
//        Factory function for creating scaled ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
//        The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
//        This function creates an instance of an `EphemerisSettings` derived `ScaledEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        unscaled_ephemeris_settings : EphemerisSettings
//            Sets base settings of ephemeris to be scaled.
//        scaling_vector : numpy.ndarray
//            Vector containing scaling factors to be applied to each element of the Cartesian state.
//        is_scaling_absolute : bool, default=false
//            Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
//
//        Returns
//        -------
//        ScaledEphemerisSettings
//
//
//    )";
//
//
//    } else if(name == "scaled" && variant==2) {
//            return R"(
//
//        Factory function for creating scaled ephemeris model settings.
//
//        Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
//        The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).
//        This function creates an instance of an `EphemerisSettings` derived `ScaledEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        unscaled_ephemeris_settings : EphemerisSettings
//            Sets base settings of ephemeris to be scaled.
//        scaling_vector_function : function < numpy.ndarray >
//            Function returning a vector with the scaling factors to be applied to each element of the Cartesian state.
//        is_scaling_absolute : bool, default=false
//            Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
//
//        Returns
//        -------
//        ScaledEphemerisSettings
//
//
//    )";
//
//
//
//    } else if(name == "tabulated" && variant==0) {
//            return R"(
//
//        Factory function for creating ephemeris model settings from tabulated data.
//
//        Factory function for settings object, defining ephemeris model to be created from tabulated data.
//        Currently the data that is provided gets interpolated by a 6th order Lagrange interpolator (hardcoded).
//        At the edges of the interpolation interval a cubic spline interpolator is used to suppres the influence of Runge's phenomenon.
//        This function creates an instance of an `EphemerisSettings` derived `TabulatedEphemerisSettings` object.
//
//
//        Parameters
//        ----------
//        body_state_history : dict
//            Dictionary of the discrete state history data from which ephemeris is to be created. Keys representing the time (float) and values representing Cartesian states (numpy.ndarray).
//        frame_origin : str, default='SSB'
//            Origin of frame in which ephemeris data is defined.
//        frame_orientation : str, default='ECLIPJ2000'
//            Orientation of frame in which ephemeris data is defined.
//
//        Returns
//        -------
//        TabulatedEphemerisSettings
//
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//
//
//
//
//namespace gravity_field {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//    } else if(name == "GravityFieldType") {
//         return R"(
//
//        Enumeration of gravity field types.
//
//        Enumeration of gravity field types supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "SphericalHarmonicsModel") {
//         return R"(
//
//        Enumeration of spherical harmonics models.
//
//        Enumeration of spherical harmonics models supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "GravityFieldSettings") {
//         return R"(
//
//        Base class for providing settings for automatic gravity field model creation.
//
//        This class is a functional base class for settings of gravity field models that require no information in addition to their type.
//        Gravity field model classes requiring additional information must be created using an object derived from this class.
//
//     )";
//
//
//    } else if(name == "GravityFieldSettings.__init__" && variant==0) {
//            return R"(
//
//    )";
//
//
//
//
//    } else if(name == "CentralGravityFieldSettings") {
//         return R"(
//
//        `GravityFieldSettings` derived class defining settings of point mass gravity field.
//
//        Derived class of `GravityFieldSettings` for central gravity fields, which are defined by a single gravitational parameter.
//
//     )";
//
//
//
//    } else if(name == "SphericalHarmonicsGravityFieldSettings") {
//         return R"(
//
//        `GravityFieldSettings` derived class defining settings of spherical harmonic gravity field representation.
//
//        Derived class of `GravityFieldSettings` for gravity fields, which are defined by a spherical harmonic gravity field representation.
//
//     )";
//
//
//
//
//    } else if(name == "central" && variant==0) {
//            return R"(
//
//        Factory function for central gravity field settings object.
//
//        Factory function for settings object, defining a point-mass gravity field model with user-defined gravitational parameter.
//        This function returns a `GravityFieldSettings` derived `CentralGravityFieldSettings` object.
//
//
//        Parameters
//        ----------
//        gravitational_parameter : float
//            None
//
//        Returns
//        -------
//        CentralGravityFieldSettings
//            `CentralGravityFieldSettings` object defined by the provided gravitational parameter.
//
//    )";
//
//
//
//    } else if(name == "central_spice" && variant==0) {
//            return R"(
//
//        Factory function to create central gravity field settings from Spice settings.
//
//        Factory function for settings object, defining a point-mass gravity field model with gravitational parameter from Spice.
//        This function returns a `GravityFieldSettings` object of gravity field type ``central_spice``.
//
//
//        Parameters
//        ----------
//        None
//            None
//
//        Returns
//        -------
//        GravityFieldSettings
//            `GravityFieldSettings` object defined by gravitational parameters from Spice settings.
//
//    )";
//
//
//
//    } else if(name == "spherical_harmonic" && variant==0) {
//            return R"(
//
//        Factory function for creating a spherical harmonics gravity field settings object.
//
//        Factory function for settings object, defining a gravity field model through spherical harmonic expansion.
//        The associated reference frame must presently be the same frame ID as the target frame of the body’s rotation model.
//        It represents the frame in which the spherical harmonic field is defined.
//        Spherical harmonic coefficients used for this environment model must *always* be fully normalized.
//        To normalize unnormalized spherical harmonic coefficients, see `spherical_harmonics_normalization`.
//        This function returns a `GravityFieldSettings` derived `SphericalHarmonicsGravityFieldSettings` object.
//
//
//        Parameters
//        ----------
//        gravitational_parameter : float
//            Gravitational parameter of gravity field.
//        reference_radius : float
//            Reference radius of spherical harmonic field expansion.
//        normalized_cosine_coefficients : numpy.ndarray
//            Cosine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.
//            As such, note that entry (0,0) of cosine coefficients should be equal to 1.
//
//        normalized_sine_coefficients : numpy.ndarray
//            Sine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.
//        associated_reference_frame : str
//            Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
//
//        Returns
//        -------
//        SphericalHarmonicsGravityFieldSettings
//            `SphericalHarmonicsGravityFieldSettings` object defined by the provided parameters.
//
//    )";
//
//
//
//    } else if(name == "spherical_harmonic_triaxial_body" && variant==0) {
//            return R"(
//
//        Factory function for spherical harmonics gravity field settings object from triaxial ellipsoid parameters.
//
//        Factory function for settings object, defining a gravity field model through spherical harmonic expansion of a homogeneous triaxial ellipsoid.
//        The constant mass distribution in the specified ellipsoid shape is expanded to obtain a spherical harmonic coefficient representation.
//        Gravity fields from this setting object are expressed in normalized spherical harmonic coefficients.
//        (%!) The x-, y- and z-axis of the ... ABC (?).
//        This function returns a `GravityFieldSettings` derived `SphericalHarmonicsGravityFieldSettings` object.
//
//
//        Parameters
//        ----------
//        axis_a : float
//            Dimension of largest axis of triaxial ellipsoid.
//        axis_b : float
//            Dimension of intermediate axis of triaxial ellipsoid.
//        axis_c : float
//            Dimension of smallest axis of triaxial ellipsoid.
//        density : float
//            Density of ellipsoid.
//        maximum_degree : int
//            Maximum degree of spherical harmonics expansion.
//        maximum_order : int
//            Maximum order of spherical harmonics expansion.
//        associated_reference_frame : str
//            Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
//
//        Returns
//        -------
//        SphericalHarmonicsGravityFieldSettings
//            `SphericalHarmonicsGravityFieldSettings` object defined by expansion of homogeneous triaxial ellipsoid.
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//
//
//
//
//namespace rotation_model {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//    } else if(name == "RotationModelType") {
//         return R"(
//
//        Enumeration of rotation model types.
//
//        Enumeration of rotation model types supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "IAUConventions") {
//         return R"(
//
//        Enumeration of IAU conventions for Earth rotation.
//
//        Enumeration of IAU conventions for Earth rotation supported by tudat.
//
//     )";
//
//
//    } else if(name == "IAUConventions.iau_2000_a" && variant==0) {
//            return R"(
//
//    )";
//
//
//
//    } else if(name == "IAUConventions.iau_2000_b" && variant==0) {
//            return R"(
//
//    )";
//
//
//
//    } else if(name == "IAUConventions.iau_2006" && variant==0) {
//            return R"(
//
//    )";
//
//
//
//
//    } else if(name == "RotationModelSettings") {
//         return R"(
//
//        Base class for providing settings for automatic rotation model creation.
//
//        This class is a functional base class for settings of rotation models that require no information in addition to their type.
//        Basic rotation model has constant orientation of the rotation axis (body-fixed z-axis) and constant rotation rate about this axis.
//        Rotation models requiring additional information must be created using the factory functions which create the specific object derived from this base class.
//
//     )";
//
//
//
//
//    } else if(name == "simple" && variant==0) {
//            return R"(
//
//        Factory function for creating simple rotation model settings.
//
//        Factory function for settings object, defining a basic rotation model with constant orientation of the rotation axis and constant rotation rate about this axis.
//        Rotation from original (inertial) to target (body-fixed) frame at some reference time ``initial_time`` is defined by the ``initial_orientation`` rotation matrix.
//        Rotation about the body-fixed z-axis is defined by the ``rotation_rate`` float variable (in rad/s).
//        This function creates an instance of a `RotationModelSettings` derived `simpleRotationModelSettings` object.
//
//
//        Parameters
//        ----------
//        base_frame : str
//            Base frame of rotation model.
//        target_frame : str
//            Target frame of rotation model.
//        initial_orientation : numpy.ndarray
//            Orientation of target frame in base frame at initial time.
//        initial_time : float
//            Initial time (reference epoch for rotation matrices).
//        rotation_rate : float
//            Constant rotation rate [rad/s] about rotational axis.
//
//        Returns
//        -------
//        SimpleRotationModelSettings
//            Simple rotation model settings object (derived from `RotationModelSettings` base class).
//
//    )";
//
//
//
//    } else if(name == "simple_from_spice" && variant==0) {
//            return R"(
//
//        Factory function for creating simple rotation model settings using initial orientation and rotaton rates from Spice.
//
//        Factory function for settings object, defining a simple rotation model with the added functionality that the initial orientation and rotation rate are extracted from Spice, as opposed to provided manually.
//        Note that *only* the initial orientation and rotation rate ( at the time defined by `initial_time` ) are extracted from Spice.
//        The distinction between the two target frame inputs is the following
//
//          - the ``target_frame`` parameter is the name of frame that Tudat assigns to the body-fixed frame
//          - the ``target_frame_spice`` is the name of the frame in Spice for which the initial orientation and rotation rate are extracted.
//
//        This function creates an instance of a `RotationModelSettings` derived `simpleRotationModelSettings` object.
//
//
//        Parameters
//        ----------
//        base_frame : str
//            Base frame of rotation model.
//        target_frame : str
//            Target frame of rotation model.
//        target_frame_spice : str
//            Spice reference of target frame.
//        initial_time : float
//            Initial time (reference epoch for rotation matrices).
//
//        Returns
//        -------
//        SimpleRotationModelSettings
//            Simple rotation model settings object (derived from RotationModelSettings base class) with target frame info inferred from Spice.
//
//    )";
//
//
//
//    } else if(name == "synchronous" && variant==0) {
//            return R"(
//
//        Factory function for creating synchronous rotational ephemeris settings.
//
//        Factory function for settings object, defining a synchronous rotation model where rotation of a body is defined from its relative orbit w.r.t. some central body. Specifically
//        - the body-fixed x-axis is *always* pointing towards the central body
//        - the body-fixed z-axis is *always* perpendicular to the orbital plane (along the direction of
//        .. math:: \mathbf{x} \cross \mathbf{v} )
//        - the body-fixed y-axis completes the right-handed reference frame
//
//        Such a model can be useful for, for instance, approximate rotation of tidally locked natural satellites or nadir-pointing spacraft.
//        This function creates an instance of a `RotationModelSettings` derived `SynchronousRotationModelSettings` object.
//
//
//        Parameters
//        ----------
//        central_body_name : str
//            Base frame of rotation model.
//        base_frame : str
//            Base frame of rotation model.
//        target_frame : str
//            Spice reference of target frame.
//
//        Returns
//        -------
//        SynchronousRotationModelSettings
//            Synchonous rotation model settings object (derived from RotationModelSettings base class).
//
//    )";
//
//
//
//    } else if(name == "spice" && variant==0) {
//            return R"(
//
//        Factory function for creating rotation model settings from the Spice interface.
//
//        Factory function for settings object, defining a rotation model directly (and entirely) from Spice interface.
//        This function creates an instance of a `RotationModelSettings` object.
//
//
//        Parameters
//        ----------
//        base_frame : str
//            Base frame of rotation model.
//        target_frame : str
//            Target frame of rotation model.
//
//        Returns
//        -------
//        RotationModelSettings
//            Rotation model settings object inferred from Spice rotational model.
//
//    )";
//
//
//
//    } else if(name == "gcrs_to_itrs" && variant==0) {
//            return R"(
//
//        Factory function for creating high-accuracy Earth rotation model settings.
//
//        Factory function for settings object, defining high-accuracy Earth rotation model according to the IERS 2010 Conventions.
//        This settings class has various options to deviate from the default settings, typical applications will use default.
//        Note that for this model the original frame must be J2000 or GCRS (in the case of the former, the frame bias between GCRS and J2000 is automatically corrected for). The target frame (e.g. body-fixed frame) name is ITRS.
//        The precession-nutation theory may be `iau_2000a` / `iau_2000b` or `iau_2006`, as implemented in the SOFA toolbox. Alternative options to modify the input (not shown above) include the EOP correction file, input time scale, short period UT1 and polar motion variations.
//        The target frame (e.g. body-fixed frame) name is ITRS.
//        This function creates an instance of a `RotationModelSettings` derived `gcrsToItrsRotationModelSettings` object.
//
//
//        Parameters
//        ----------
//        precession_nutation_theory : default=tba::iau_2006
//            Setting theory for modelling Earth nutation.
//
//        base_frame : str, default='GCRS'
//            Base frame of rotation model
//
//        Returns
//        -------
//        GcrsToItrsRotationModelSettings
//            High-accuracy Earth rotation model settings object (derived from RotationModelSettings base class).
//
//    )";
//
//
//
//    } else if(name == "constant" && variant==0) {
//            return R"(
//
//        Factory function for creating simple rotation model settings for target-frames with constant orientation.
//
//        Factory function for settings object, defining simple rotation model setting objects with constant rotation matrix.
//        These model settings are for target frames which do not have a rotational rate in the base frame and are fully defined by their initial orientation.
//        This function creates an instance of a `RotationModelSettings` derived `SimpleRotationModelSettings` object.
//
//
//        Parameters
//        ----------
//        base_frame : str
//            Base frame of rotation model.
//        target_frame : str
//            Target frame of rotation model.
//        initial_orientation : numpy.ndarray
//            Orientation of target frame in base frame at initial time (constant throughout).
//
//        Returns
//        -------
//        SimpleRotationModelSettings
//            Simple rotation model settings object (derived from RotationModelSettings base class) with constant orientation of target in base frame.
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//
//
//
//}
//
//
//
//
//
//namespace propagation_setup {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//
//namespace acceleration {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//    } else if(name == "AvailableAcceleration") {
//         return R"(
//
//        Enumeration of available acceleration types.
//
//        Enumeration of acceleration types supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "ThrustFrames") {
//         return R"(
//
//        Enumeration of available thrust frame types.
//
//        Enumeration of thrust frame types supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "ThrustMagnitudeTypes") {
//         return R"(
//
//        Enumeration of available thrust magnitude types.
//
//        Enumeration of thrust magnitude types supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "ThrustDirectionTypes") {
//         return R"(
//
//        Enumeration of available thrust direction types.
//
//        Enumeration of thrust direction types supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "AccelerationSettings") {
//         return R"(
//
//        Functional base class to define settings for accelerations.
//
//     )";
//
//
//
//    } else if(name == "SphericalHarmonicAccelerationSettings") {
//         return R"(
//
//        `AccelerationSettings`-derived class to define settings for the spherical harmonic acceleration.
//
//     )";
//
//
//
//    } else if(name == "MutualSphericalHarmonicAccelerationSettings") {
//         return R"(
//
//        `AccelerationSettings`-derived class to define settings for the mutual spherical harmonic acceleration.
//
//     )";
//
//
//
//    } else if(name == "RelativisticAccelerationCorrectionSettings") {
//         return R"(
//
//        `AccelerationSettings`-derived class to define settings for the relativistic acceleration correction.
//
//     )";
//
//
//
//    } else if(name == "EmpiricalAccelerationSettings") {
//         return R"(
//
//        `AccelerationSettings`-derived class to define settings for the empirical acceleration.
//
//     )";
//
//
//
//    } else if(name == "CustomAccelerationSettings") {
//         return R"(
//
//        `AccelerationSettings`-derived class to define settings for custom acceleration.
//
//     )";
//
//
//
//    } else if(name == "DirectTidalDissipationAccelerationSettings") {
//         return R"(
//
//        `AccelerationSettings`-derived class to define settings for direct tidal dissipation acceleration.
//
//     )";
//
//
//
//    } else if(name == "ThrustAccelerationSettings") {
//         return R"(
//
//        `AccelerationSettings`-derived class to define settings for direct tidal dissipation acceleration.
//
//     )";
//
//
//
//    } else if(name == "MomentumWheelDesaturationAccelerationSettings") {
//         return R"(
//
//        `AccelerationSettings`-derived class to define settings for momentum wheel desaturation acceleration.
//
//     )";
//
//
//
//    } else if(name == "ThrustDirectionSettings") {
//         return R"(
//
//        Functional base class to define settings for the thrust direction.
//
//     )";
//
//
//
//    } else if(name == "ThrustDirectionFromStateGuidanceSettings") {
//         return R"(
//
//        `ThrustDirectionSettings`-derived class to define settings for the thrust direction from the current state.
//
//     )";
//
//
//
//    } else if(name == "CustomThrustDirectionSettings") {
//         return R"(
//
//        `ThrustDirectionSettings`-derived class to define settings for a custom thrust direction.
//
//     )";
//
//
//
//    } else if(name == "CustomThrustOrientationSettings") {
//         return R"(
//
//        `ThrustDirectionSettings`-derived class to define settings for a custom thrust orientation.
//
//     )";
//
//
//
//    } else if(name == "MeeCostateBasedThrustDirectionSettings") {
//         return R"(
//
//        `ThrustDirectionSettings`-derived class to define settings for the thrust direction from Modified Equinoctial Elements (MEE) costates.
//
//     )";
//
//
//
//    } else if(name == "ThrustMagnitudeSettings") {
//         return R"(
//
//        Functional base class to define settings for the thrust magnitude.
//
//     )";
//
//
//
//    } else if(name == "ConstantThrustMagnitudeSettings") {
//         return R"(
//
//        `ThrustMagnitudeSettings`-derived class to define settings for constant thrust magnitude.
//
//     )";
//
//
//
//    } else if(name == "FromFunctionThrustMagnitudeSettings") {
//         return R"(
//
//        `ThrustMagnitudeSettings`-derived class to define settings for constant thrust magnitude.
//
//     )";
//
//
//
//
//    } else if(name == "point_mass_gravity" && variant==0) {
//            return R"(
//
//        Creates settings for the point-mass gravity acceleration.
//
//        Creates settings for the point-mass gravity acceleration. The body exerting the acceleration needs to have a
//        gravity field model defined.
//
//
//        Returns
//        -------
//        AccelerationSettings
//            Acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "aerodynamic" && variant==0) {
//            return R"(
//
//        Creates settings for the aerodynamic acceleration.
//
//        Creates settings for the aerodynamic acceleration. The body exerting the acceleration needs to have an
//        atmosphere defined.
//
//
//        Returns
//        -------
//        AccelerationSettings
//            Acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "cannonball_radiation_pressure" && variant==0) {
//            return R"(
//
//        Creates settings for the cannonball radiation pressure acceleration.
//
//        Creates settings for the radiation pressure acceleration, for which a cannonball model is used. In this model,
//        the effective acceleration is colinear with the vector connecting the source of radiation and the target.
//        The body undergoing the acceleration needs to have a radiation pressure model defined, while the body emitting
//        radiation needs to have radiative properties defined (the Sun has default ones).
//
//
//        Returns
//        -------
//        AccelerationSettings
//            Acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "spherical_harmonic_gravity" && variant==0) {
//            return R"(
//
//        Creates settings for the spherical harmonic gravity acceleration.
//
//        Creates settings for the spherical harmonic gravity acceleration, accounting for a finite (given) number
//        of degree and order. The body exerting the acceleration needs to have a spherical harmonic gravity field model
//        defined.
//
//
//        Parameters
//        ----------
//        maximum_degree : int
//            Maximum degree of the spherical harmonic expansion.
//        maximum_order : int
//            Maximum order of the spherical harmonic expansion.
//
//        Returns
//        -------
//        SphericalHarmonicAccelerationSettings
//            Spherical harmonic acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "mutual_spherical_harmonic_gravity" && variant==0) {
//            return R"(
//
//        Creates settings for the mutual spherical harmonic gravity acceleration.
//
//        Creates settings for the mutual spherical harmonic gravity acceleration, accounting for a finite (given) number
//        of degree and order for both bodies. Both the body exerting the acceleration and the body undergoing it need to
//        have spherical harmonic gravity field models defined. In addition, the body undergoing the acceleration needs to
//        have a rotational model defined. For the case where a third-body mutual spherical harmonic acceleration,
//        additional parameters have to be provided that denote the expansion degree/order of the central body.
//
//
//        Parameters
//        ----------
//        maximum_degree_body_exerting : int
//            Maximum degree of the spherical harmonic expansion for the body exerting the acceleration.
//        maximum_order_body_exerting : int
//            Maximum order of the spherical harmonic expansion for the body exerting the acceleration.
//        maximum_degree_body_undergoing : int
//            Maximum degree of the spherical harmonic expansion for the body undergoing the acceleration.
//        maximum_order_body_undergoing : int
//            Maximum order of the spherical harmonic expansion for the body undergoing the acceleration.
//        maximum_degree_central_body : int, default=0
//            Maximum degree of the spherical harmonic expansion for the central body, if needed.
//        maximum_order_central_body : int, default=0
//            Maximum order of the spherical harmonic expansion for the central body, if needed.
//
//        Returns
//        -------
//        MutualSphericalHarmonicAccelerationSettings
//            Spherical harmonic acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "relativistic_correction" && variant==0) {
//            return R"(
//
//        Creates settings for the relativistic acceleration correction.
//
//        Creates settings for typical relativistic acceleration corrections: the Schwarzschild, Lense-Thirring and de
//        Sitter terms (see 'General relativity and Space Geodesy' by L. Combrinck, 2012). It implements the model of
//        2010 Conventions (chapter 10, section 3). Here, the ‘primary body’ for a planetary orbiter should always be set
//        as the Sun (only relevant for de Sitter correction). The angular momentum vector of the orbited body is only
//        relevant for Lense-Thirring correction.
//
//
//        Parameters
//        ----------
//        use_schwarzschild : bool
//            Maximum degree of the spherical harmonic expansion for the body exerting the acceleration.
//        use_lense_thirring : bool
//            Maximum order of the spherical harmonic expansion for the body exerting the acceleration.
//        use_de_sitter : bool
//            Maximum degree of the spherical harmonic expansion for the body undergoing the acceleration.
//        de_sitter_central_body : str, default=""
//            Maximum order of the spherical harmonic expansion for the body undergoing the acceleration.
//        lense_thirring_angular_momentum : numpy.ndarray, default=numpy.array([0, 0, 0])
//            Maximum degree of the spherical harmonic expansion for the central body, if needed.
//
//        Returns
//        -------
//        RelativisticAccelerationCorrectionSettings
//            Relativistic acceleration correction settings object.
//
//    )";
//
//
//
//    } else if(name == "empirical" && variant==0) {
//            return R"(
//
//        Creates settings for empirical acceleration.
//
//        Creates settings for empirical accelerations. These are expressed in the
//        RSW frame, for which the mangitude is determined empirically (typically during an orbit determination process).
//        The acceleration components are defined according to Montenbruck and Gill (2000), with a total of 9 components:
//        a constant, sine and cosine term (with true anomaly as argument) for each of the three independent directions of
//        the RSW frame.
//
//
//        Parameters
//        ----------
//        constant_acceleration : numpy.ndarray, default=numpy.array([0, 0, 0])
//            Constant term, defined in the RSW frame.
//        sine_acceleration : numpy.ndarray, default=numpy.array([0, 0, 0])
//            Sine term (function of the true anomaly), defined in the RSW frame..
//        cosine_acceleration : numpy.ndarray, default=numpy.array([0, 0, 0])
//            Cosine term (function of the true anomaly), defined in the RSW frame..
//
//        Returns
//        -------
//        EmpiricalAccelerationSettings
//            Empirical acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "custom" && variant==0) {
//            return R"(
//
//        Creates settings for custom acceleration.
//
//        Creates settings for empirical accelerations. These are expressed in the
//        RSW frame, for which the mangitude is determined empirically (typically during an orbit determination process).
//        The acceleration components are defined according to Montenbruck and Gill (2000), with a total of 9 components:
//        a constant, sine and cosine term (with true anomaly as argument) for each of the three independent directions of
//        the RSW frame.
//
//
//        Parameters
//        ----------
//        acceleration_function : Callable[[float], list]
//            Custom acceleration function with time as an independent variable.
//        scaling_function : Callable[[float], float], default=None
//            Scaling function with time as an independent variable to be multiplied by the custom acceleration function.
//
//        Returns
//        -------
//        CustomAccelerationSettings
//            Custom acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "direct_tidal_dissipation_acceleration" && variant==0) {
//            return R"(
//
//        Creates settings for custom acceleration.
//
//        Creates settings for tidal accelerations. The direct of tidal effects in a satellite system is applied directly as
//        an acceleration (as opposed to a modification of spherical harmonic coefficients).
//        The model is based on Lainey et al. (2007, 2012). It can compute the acceleration due to tides, and in
//        particular tidal dissipation, on a planetary satellite. The acceleration computed can account for either the
//        effect of tide raised on the satellite by the planet or on the planet by the satellite. The satellite is assumed
//        to be tidally locked to the planet.
//
//
//        Parameters
//        ----------
//        k2_love_number : float
//            Value of the k2 Love number.
//        time_lag : float
//            Value of the tidal time lag.
//        include_direct_radial_component : bool, default=True
//            It denotes whether the term independent of the time lag is to be computed.
//        use_tide_raised_on_planet : bool, default=True
//            It denotes whether the tide raised on the planet is to be modelled (if true) or the tide raised on the satellite (if false).
//
//        Returns
//        -------
//        DirectTidalDissipationAccelerationSettings
//            Direct tidal dissipation acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "momentum_wheel_desaturation_acceleration" && variant==0) {
//            return R"(
//
//        Creates settings for momentum wheel desaturation acceleration.
//
//        The acceleration model is purpose-built to represent short bursts of thrust, such as a momentum wheel desaturation.
//        A typical use case is precise orbit determination, but the functionality can be used just as well in propagation
//        (for instance to model an impulsive manuever in a continuous manner when going from preliminary modelling to
//        'full' modelling). The thrust is modelled similarly to Fig. 3 of Alessi et al. (2012), with the main difference
//        being that a third-order polynomial to go from zero acceleration to the maximum acceleration level is employed.
//        By using a 3rd-order polynomial and imposing continuity in the value and first derivative of the acceleration,
//        defining the 'rise time' (time it takes acceleration to go from 0 to its maximum level), the total time where
//        there is non-zero thrust ('total maneuver time'), and the total Delta V exerted by a single maneuver,
//        the acceleration profile is fully defined.
//
//
//        Parameters
//        ----------
//        thrust_mid_times : list[float]
//            Set of middle point in times in the maneuver denoting the epoch of each maneuver.
//        delta_v_values : list[numpy.ndarray]
//            Set of delta V, one for each maneuver.
//        total_maneuver_time : float
//            Total duration of every maneuver.
//        maneuver_rise_time : float
//            Time taken by the acceleration to go from zero to its maximum level.
//
//        Returns
//        -------
//        MomentumWheelDesaturationAccelerationSettings
//            Momentum wheel desaturation acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "thrust_acceleration" && variant==0) {
//            return R"(
//
//        Creates settings for thrust acceleration from thrust guidance settings.
//
//        Creates settings for thrust acceleration from thrust guidance settings. The thrust direction and magnitude are
//        supplied in the form of dedicated settings objects (see the API for the respective classes).
//
//
//        Parameters
//        ----------
//        thrust_direction_settings : ThrustDirectionSettings
//            Thrust direction settings object.
//        thrust_magnitude_settings : ThrustMagnitudeSettings
//            Thrust magnitude settings object.
//
//        Returns
//        -------
//        ThrustAccelerationSettings
//            Thrust acceleration settings object.
//
//    )";
//
//
//    } else if(name == "thrust_acceleration" && variant==1) {
//            return R"(
//
//        Creates settings for thrust acceleration from interpolated thrust data with variable magnitude.
//
//        Creates settings for thrust acceleration from interpolated thrust data. The thrust direction is defined through
//        the related interpolator (which uses time as independent variable) and it returns the thrust direction vector in
//        the specified frame (it can be local or inertial). The variable thrust magnitude is computed from the specific impulse, given as a function of time.
//
//
//        Parameters
//        ----------
//        data_interpolation_settings : DataInterpolationSettings<float, numpy.ndarray>
//            Interpolator object that provides the thrust direction vector in the given thrust frame as a function of time.
//        specific_impulse_function : Callable[[double], double]
//            Specific impulse provided as a function of time.
//        thrust_frame : ThrustFrames, default=unspecified_thrust_frame
//            Frame in which the thrust direction vector is represented.
//        central_body : str, default=""
//            Central body that is the origin of the thrust frame (if different from the vehicle, otherwise empty by default).
//
//        Returns
//        -------
//        ThrustAccelerationSettings
//            Thrust acceleration settings object.
//
//    )";
//
//
//    } else if(name == "thrust_acceleration" && variant==2) {
//            return R"(
//
//        Creates settings for thrust acceleration from interpolated thrust data with constant magnitudee.
//
//        Creates settings for thrust acceleration from interpolated thrust data. The thrust direction is defined through
//        the related interpolator (which uses time as independent variable) and it returns the thrust direction vector in
//        the specified frame (it can be local or inertial). The constant thrust magnitude is computed from the constant
//        specific impulse.
//
//
//        Parameters
//        ----------
//        data_interpolation_settings : DataInterpolationSettings<float, numpy.ndarray>
//            Interpolator object that provides the thrust direction vector in the given thrust frame as a function of time.
//        constant_specific_impulse : float
//            Constant specific impulse.
//        thrust_frame : ThrustFrames, default=unspecified_thrust_frame
//            Frame in which the thrust direction vector is represented.
//        central_body : str, default=""
//            Central body that is the origin of the thrust frame (if different from the vehicle, otherwise empty by default).
//
//        Returns
//        -------
//        ThrustAccelerationSettings
//            Thrust acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "thrust_acceleration" && variant==0) {
//            return R"(
//
//        Creates settings for thrust acceleration from thrust guidance settings.
//
//        Creates settings for thrust acceleration from thrust guidance settings. The thrust direction and magnitude are
//        supplied in the form of dedicated settings objects (see the API for the respective classes).
//
//
//        Parameters
//        ----------
//        thrust_direction_settings : ThrustDirectionSettings
//            Thrust direction settings object.
//        thrust_magnitude_settings : ThrustMagnitudeSettings
//            Thrust magnitude settings object.
//
//        Returns
//        -------
//        ThrustAccelerationSettings
//            Thrust acceleration settings object.
//
//    )";
//
//
//    } else if(name == "thrust_acceleration" && variant==1) {
//            return R"(
//
//        Creates settings for thrust acceleration from interpolated thrust data with variable magnitude.
//
//        Creates settings for thrust acceleration from interpolated thrust data. The thrust direction is defined through
//        the related interpolator (which uses time as independent variable) and it returns the thrust direction vector in
//        the specified frame (it can be local or inertial). The variable thrust magnitude is computed from the specific impulse, given as a function of time.
//
//
//        Parameters
//        ----------
//        data_interpolation_settings : DataInterpolationSettings<float, numpy.ndarray>
//            Interpolator object that provides the thrust direction vector in the given thrust frame as a function of time.
//        specific_impulse_function : Callable[[double], double]
//            Specific impulse provided as a function of time.
//        thrust_frame : ThrustFrames, default=unspecified_thrust_frame
//            Frame in which the thrust direction vector is represented.
//        central_body : str, default=""
//            Central body that is the origin of the thrust frame (if different from the vehicle, otherwise empty by default).
//
//        Returns
//        -------
//        ThrustAccelerationSettings
//            Thrust acceleration settings object.
//
//    )";
//
//
//    } else if(name == "thrust_acceleration" && variant==2) {
//            return R"(
//
//        Creates settings for thrust acceleration from interpolated thrust data with constant magnitudee.
//
//        Creates settings for thrust acceleration from interpolated thrust data. The thrust direction is defined through
//        the related interpolator (which uses time as independent variable) and it returns the thrust direction vector in
//        the specified frame (it can be local or inertial). The constant thrust magnitude is computed from the constant
//        specific impulse.
//
//
//        Parameters
//        ----------
//        data_interpolation_settings : DataInterpolationSettings<float, numpy.ndarray>
//            Interpolator object that provides the thrust direction vector in the given thrust frame as a function of time.
//        constant_specific_impulse : float
//            Constant specific impulse.
//        thrust_frame : ThrustFrames, default=unspecified_thrust_frame
//            Frame in which the thrust direction vector is represented.
//        central_body : str, default=""
//            Central body that is the origin of the thrust frame (if different from the vehicle, otherwise empty by default).
//
//        Returns
//        -------
//        ThrustAccelerationSettings
//            Thrust acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "thrust_acceleration" && variant==0) {
//            return R"(
//
//        Creates settings for thrust acceleration from thrust guidance settings.
//
//        Creates settings for thrust acceleration from thrust guidance settings. The thrust direction and magnitude are
//        supplied in the form of dedicated settings objects (see the API for the respective classes).
//
//
//        Parameters
//        ----------
//        thrust_direction_settings : ThrustDirectionSettings
//            Thrust direction settings object.
//        thrust_magnitude_settings : ThrustMagnitudeSettings
//            Thrust magnitude settings object.
//
//        Returns
//        -------
//        ThrustAccelerationSettings
//            Thrust acceleration settings object.
//
//    )";
//
//
//    } else if(name == "thrust_acceleration" && variant==1) {
//            return R"(
//
//        Creates settings for thrust acceleration from interpolated thrust data with variable magnitude.
//
//        Creates settings for thrust acceleration from interpolated thrust data. The thrust direction is defined through
//        the related interpolator (which uses time as independent variable) and it returns the thrust direction vector in
//        the specified frame (it can be local or inertial). The variable thrust magnitude is computed from the specific impulse, given as a function of time.
//
//
//        Parameters
//        ----------
//        data_interpolation_settings : DataInterpolationSettings<float, numpy.ndarray>
//            Interpolator object that provides the thrust direction vector in the given thrust frame as a function of time.
//        specific_impulse_function : Callable[[double], double]
//            Specific impulse provided as a function of time.
//        thrust_frame : ThrustFrames, default=unspecified_thrust_frame
//            Frame in which the thrust direction vector is represented.
//        central_body : str, default=""
//            Central body that is the origin of the thrust frame (if different from the vehicle, otherwise empty by default).
//
//        Returns
//        -------
//        ThrustAccelerationSettings
//            Thrust acceleration settings object.
//
//    )";
//
//
//    } else if(name == "thrust_acceleration" && variant==2) {
//            return R"(
//
//        Creates settings for thrust acceleration from interpolated thrust data with constant magnitudee.
//
//        Creates settings for thrust acceleration from interpolated thrust data. The thrust direction is defined through
//        the related interpolator (which uses time as independent variable) and it returns the thrust direction vector in
//        the specified frame (it can be local or inertial). The constant thrust magnitude is computed from the constant
//        specific impulse.
//
//
//        Parameters
//        ----------
//        data_interpolation_settings : DataInterpolationSettings<float, numpy.ndarray>
//            Interpolator object that provides the thrust direction vector in the given thrust frame as a function of time.
//        constant_specific_impulse : float
//            Constant specific impulse.
//        thrust_frame : ThrustFrames, default=unspecified_thrust_frame
//            Frame in which the thrust direction vector is represented.
//        central_body : str, default=""
//            Central body that is the origin of the thrust frame (if different from the vehicle, otherwise empty by default).
//
//        Returns
//        -------
//        ThrustAccelerationSettings
//            Thrust acceleration settings object.
//
//    )";
//
//
//
//    } else if(name == "get_propulsion_input_variables" && variant==0) {
//            return R"(
//
//        Function to create a list of functions that compute and return independent variables for the thrust.
//
//        Function to create a list of functions that compute and return independent variables for thrust and/or specific
//        impulse. This parameterization is used to create a specific thrust magnitude type (see thrust magnitude from
//        dependent variables). This function retrieves all input functions from the environment and a list of user-defined
//        functions.
//
//
//        Parameters
//        ----------
//        body_with_guidance : Body
//            Body object whose thrust guidance should be defined.
//        independent_variables : list[ThrustIndependentVariables]
//            Set of dependent variables that should be used to compute the thrust.
//        guidance_input_functions : list[Callable[[], float], default=[]
//            Set of functions to compute the thrust, each associated to a specific dependent variable.
//    )";
//
//
//
//    } else if(name == "thrust_direction_from_state_guidance" && variant==0) {
//            return R"(
//
//        Create thrust direction settings from the state guidance.
//
//        Factory function that creates thrust direction settings from the state guidance. In various simplified cases,
//        the thrust direction can be assumed to be in line with either the position or velocity of the body of interest
//        with respect to some other body.
//
//
//        Parameters
//        ----------
//        central_body : Body
//            Central body with respect to which the position and velocity of the body undergoing the thrust acceleration are computed.
//        is_colinear_with_velocity : bool
//            Whether the thrust direction is colinear with the velocity (true) or the position vector with respect to some other body (false).
//        direction_is_opposite_to_vector : bool
//            Whether the thrust is pointing towards the thrusting body (true) or the central body (false).
//
//        Returns
//        -------
//        ThrustDirectionFromStateGuidanceSettings
//            Thrust direction from state guidance settings object.
//
//    )";
//
//
//
//    } else if(name == "thrust_from_existing_body_orientation" && variant==0) {
//            return R"(
//
//        Create thrust direction settings from the existing body orientation.
//
//        Factory function that creates thrust direction settings from the existing body orientation. In some cases,
//        the vehicle’s orientation may be predetermined, either due to aerodynamic guidance or to the concurrent
//        propagation of the rotational equations of motion. In such a case, the thrust direction is computed from the
//        body-fixed thrust direction (defined in ThrustMagnitudeSettings) and the existing vehicle orientation.
//
//    )";
//
//
//
//    } else if(name == "custom_thrust_orientation" && variant==0) {
//            return R"(
//
//        Create custom thrust orientation settings, expressed as a rotation matrix.
//
//        Factory function that creates custom thrust orientation settings, expressed through a rotation matrix.
//        As an alternative expression for generalized thrust direction guidance, the thrust orientation can be defined as
//        an arbitrary function of time. As with the custom thrust direction, this allows a broad range of options to be
//        defined, at the expense of increased complexity (somehow the thrust orientation needs to be manually defined).
//        The thrust orientation is provided through a rotation matrix representing the rotation
//        from body-fixed thrust direction to the inertial thrust direction.
//
//
//        Parameters
//        ----------
//        thrust_orientation_function : Callable[[float], numpy.ndarray]
//            Function of time returning the matrix representing the rotation between the thrust direction in the body-fixed frame to the inertial frame.
//
//        Returns
//        -------
//        CustomThrustOrientationSettings
//            Custom thrust orientation settings object.
//
//    )";
//
//
//
//    } else if(name == "custom_thrust_direction" && variant==0) {
//            return R"(
//
//        Create custom thrust direction settings, expressed as a vector in the inertial frame.
//
//        Factory function that creates custom thrust direction settings, expressed as a unit vector in the inertial frame.
//        For a generalized thrust direction guidance, the thrust can be defined as an arbitrary function of time.
//        This allows a broad range of options to be defined, at the expense of increased complexity (somehow the thrust
//        direction needs to be manually defined).
//
//
//        Parameters
//        ----------
//        thrust_direction_function : Callable[[float], numpy.ndarray]
//            Function of time returning the thrust direction in the inertial frame.
//
//        Returns
//        -------
//        CustomThrustDirectionSettings
//            Custom thrust direction settings object.
//
//    )";
//
//
//    } else if(name == "custom_thrust_direction" && variant==1) {
//            return R"(
//
//        Create thrust direction settings, expressed through modified equinoctial elements costates.
//
//        Factory function that creates thrust direction settings, expressed through modified equinoctial elements costates.
//        By using these settings for the thrust direction, the so-called co-states of the Modified Equinoctial elements
//        are used to determine the direction of the thrust. Details of this model are given by Kluever (2010),
//        Boudestijn (2014) and Hogervorst (2017). This function takes variable costates as an interpolator over time.
//
//
//        Parameters
//        ----------
//        vehicle_name : str
//            Name of the body undergoing thrust.
//        central_body_name : str
//            Name of the central body with respect to which the Modified Equinoctial Elements are computed.
//        costate_interpolator : OneDimensionalInterpolator<float, numpy.ndarray>
//            Interpolator object returning the five costates with time as an independent variable.
//
//        Returns
//        -------
//        MeeCostateBasedThrustDirectionSettings
//            Modified Equinoctial Elements costate-based thrust direction settings object.
//
//    )";
//
//
//    } else if(name == "custom_thrust_direction" && variant==2) {
//            return R"(
//
//        Create thrust direction settings, expressed through modified equinoctial elements costates.
//
//        Factory function that creates thrust direction settings, expressed through modified equinoctial elements costates.
//        By using these settings for the thrust direction, the so-called co-states of the Modified Equinoctial elements
//        are used to determine the direction of the thrust. Details of this model are given by Kluever (2010),
//        Boudestijn (2014) and Hogervorst (2017). This function takes constant costates.
//
//
//        Parameters
//        ----------
//        vehicle_name : str
//            Name of the body undergoing thrust.
//        central_body_name : str
//            Name of the central body with respect to which the Modified Equinoctial Elements are computed.
//        constant_costates : numpy.ndarray
//            Set of five constant costates.
//
//        Returns
//        -------
//        MeeCostateBasedThrustDirectionSettings
//            Modified Equinoctial Elements costate-based thrust direction settings object.
//
//    )";
//
//
//
//    } else if(name == "custom_thrust_direction" && variant==0) {
//            return R"(
//
//        Create custom thrust direction settings, expressed as a vector in the inertial frame.
//
//        Factory function that creates custom thrust direction settings, expressed as a unit vector in the inertial frame.
//        For a generalized thrust direction guidance, the thrust can be defined as an arbitrary function of time.
//        This allows a broad range of options to be defined, at the expense of increased complexity (somehow the thrust
//        direction needs to be manually defined).
//
//
//        Parameters
//        ----------
//        thrust_direction_function : Callable[[float], numpy.ndarray]
//            Function of time returning the thrust direction in the inertial frame.
//
//        Returns
//        -------
//        CustomThrustDirectionSettings
//            Custom thrust direction settings object.
//
//    )";
//
//
//    } else if(name == "custom_thrust_direction" && variant==1) {
//            return R"(
//
//        Create thrust direction settings, expressed through modified equinoctial elements costates.
//
//        Factory function that creates thrust direction settings, expressed through modified equinoctial elements costates.
//        By using these settings for the thrust direction, the so-called co-states of the Modified Equinoctial elements
//        are used to determine the direction of the thrust. Details of this model are given by Kluever (2010),
//        Boudestijn (2014) and Hogervorst (2017). This function takes variable costates as an interpolator over time.
//
//
//        Parameters
//        ----------
//        vehicle_name : str
//            Name of the body undergoing thrust.
//        central_body_name : str
//            Name of the central body with respect to which the Modified Equinoctial Elements are computed.
//        costate_interpolator : OneDimensionalInterpolator<float, numpy.ndarray>
//            Interpolator object returning the five costates with time as an independent variable.
//
//        Returns
//        -------
//        MeeCostateBasedThrustDirectionSettings
//            Modified Equinoctial Elements costate-based thrust direction settings object.
//
//    )";
//
//
//    } else if(name == "custom_thrust_direction" && variant==2) {
//            return R"(
//
//        Create thrust direction settings, expressed through modified equinoctial elements costates.
//
//        Factory function that creates thrust direction settings, expressed through modified equinoctial elements costates.
//        By using these settings for the thrust direction, the so-called co-states of the Modified Equinoctial elements
//        are used to determine the direction of the thrust. Details of this model are given by Kluever (2010),
//        Boudestijn (2014) and Hogervorst (2017). This function takes constant costates.
//
//
//        Parameters
//        ----------
//        vehicle_name : str
//            Name of the body undergoing thrust.
//        central_body_name : str
//            Name of the central body with respect to which the Modified Equinoctial Elements are computed.
//        constant_costates : numpy.ndarray
//            Set of five constant costates.
//
//        Returns
//        -------
//        MeeCostateBasedThrustDirectionSettings
//            Modified Equinoctial Elements costate-based thrust direction settings object.
//
//    )";
//
//
//
//    } else if(name == "custom_thrust_direction" && variant==0) {
//            return R"(
//
//        Create custom thrust direction settings, expressed as a vector in the inertial frame.
//
//        Factory function that creates custom thrust direction settings, expressed as a unit vector in the inertial frame.
//        For a generalized thrust direction guidance, the thrust can be defined as an arbitrary function of time.
//        This allows a broad range of options to be defined, at the expense of increased complexity (somehow the thrust
//        direction needs to be manually defined).
//
//
//        Parameters
//        ----------
//        thrust_direction_function : Callable[[float], numpy.ndarray]
//            Function of time returning the thrust direction in the inertial frame.
//
//        Returns
//        -------
//        CustomThrustDirectionSettings
//            Custom thrust direction settings object.
//
//    )";
//
//
//    } else if(name == "custom_thrust_direction" && variant==1) {
//            return R"(
//
//        Create thrust direction settings, expressed through modified equinoctial elements costates.
//
//        Factory function that creates thrust direction settings, expressed through modified equinoctial elements costates.
//        By using these settings for the thrust direction, the so-called co-states of the Modified Equinoctial elements
//        are used to determine the direction of the thrust. Details of this model are given by Kluever (2010),
//        Boudestijn (2014) and Hogervorst (2017). This function takes variable costates as an interpolator over time.
//
//
//        Parameters
//        ----------
//        vehicle_name : str
//            Name of the body undergoing thrust.
//        central_body_name : str
//            Name of the central body with respect to which the Modified Equinoctial Elements are computed.
//        costate_interpolator : OneDimensionalInterpolator<float, numpy.ndarray>
//            Interpolator object returning the five costates with time as an independent variable.
//
//        Returns
//        -------
//        MeeCostateBasedThrustDirectionSettings
//            Modified Equinoctial Elements costate-based thrust direction settings object.
//
//    )";
//
//
//    } else if(name == "custom_thrust_direction" && variant==2) {
//            return R"(
//
//        Create thrust direction settings, expressed through modified equinoctial elements costates.
//
//        Factory function that creates thrust direction settings, expressed through modified equinoctial elements costates.
//        By using these settings for the thrust direction, the so-called co-states of the Modified Equinoctial elements
//        are used to determine the direction of the thrust. Details of this model are given by Kluever (2010),
//        Boudestijn (2014) and Hogervorst (2017). This function takes constant costates.
//
//
//        Parameters
//        ----------
//        vehicle_name : str
//            Name of the body undergoing thrust.
//        central_body_name : str
//            Name of the central body with respect to which the Modified Equinoctial Elements are computed.
//        constant_costates : numpy.ndarray
//            Set of five constant costates.
//
//        Returns
//        -------
//        MeeCostateBasedThrustDirectionSettings
//            Modified Equinoctial Elements costate-based thrust direction settings object.
//
//    )";
//
//
//
//    } else if(name == "custom_thrust_magnitude" && variant==0) {
//            return R"(
//
//        Create thrust magnitude settings from a custom thrust magnitude function.
//
//        Factory function that creates constant thrust magnitude settings. The specific impulse to use for the thrust is
//        also supplied when applying a mass rate model in the propagation of the vehicle dynamics, relating the thrust
//        to the mass decrease of the vehicle.
//
//
//        Parameters
//        ----------
//        thrust_magnitude : float
//            Value of the constant thrust magnitude.
//        specific_impulse : float
//            Value of the constant specific impulse, used to link the thrust model to the mass propagation.
//        body_fixed_thrust_direction : numpy.ndarray, default=numpy.ndarray([])
//            Constant body-fixed thrust direction (positive x-direction by default). Note that this should be a unit-vector representing the direction opposite to the nozzle direction.
//
//        Returns
//        -------
//        ConstantThrustMagnitudeSettings
//            Constant thrust magnitude settings object.
//
//    )";
//
//
//    } else if(name == "custom_thrust_magnitude" && variant==1) {
//            return R"(
//
//        Create thrust magnitude settings from a custom thrust magnitude function.
//
//        Factory function that creates thrust magnitude from a custom thrust magnitude function.
//        This model defines a thrust force and specific impulse that can vary with time. The specific impulse is also
//        provided to apply a mass rate model in the propagation the vehicle dynamics, relating the thrust to the mass
//        decrease of the vehicle. Note that, if you wish to use a constant value for any or all of the first three
//        arguments, lambda expression can be used. Presently, the definition of the thrust direction in the body-fixed
//        frame is also defined through these derived classes. In essence, the ThrustMagnitudeSettings defines all local
//        (to the vehicle systems) settings for the thrust, while ThrustDirectionGuidanceSettings defines how the full
//        vehicle must orient itself in space for the required thrust direction to be achieved. At present, there is no
//        direct option for thrust-vector control (i.e. modifying the thrust direction in the body-fixed frame).
//
//
//        Parameters
//        ----------
//        thrust_magnitude_function : Callable[[float], float]
//            Function of time returning the value of the thrust magnitude.
//        specific_impulse_function : Callable[[float], float]
//            Function of time returning the value of the specific impulse, useful to link the mass propagation to the thrust model.
//        is_engine_on_function : Callable[[float], bool], default=lambda t: true
//            Function of time returning a boolean, denoting  whether the thrust should be engaged at all (e.g. thrust is 0 N if it returns false). It is useful to link the mass propagation to the thrust model.
//        body_fixed_thrust_direction
//            None
//        Callable[[], numpy.ndarray], default=lambda t: numpy.ndarray([])
//            Constant body-fixed thrust direction (positive x-direction by default). Note that this function should be a unit-vector representing the direction opposite to the nozzle direction. This setting can be used to incorporate thrust-vector control (TVC) into the thrust.
//        custom_thrust_reset_function : Callable[[float], ], default=lambda t: None
//            Function of time that updates any relevant aspects of the environment/system models, called before retrieving the thrust magnitude, specific impulse, and body-fixed thrust direction.
//
//        Returns
//        -------
//        FromFunctionThrustMagnitudeSettings
//            From function thrust magnitude settings object.
//
//    )";
//
//
//
//    } else if(name == "custom_thrust_magnitude" && variant==0) {
//            return R"(
//
//        Create thrust magnitude settings from a custom thrust magnitude function.
//
//        Factory function that creates constant thrust magnitude settings. The specific impulse to use for the thrust is
//        also supplied when applying a mass rate model in the propagation of the vehicle dynamics, relating the thrust
//        to the mass decrease of the vehicle.
//
//
//        Parameters
//        ----------
//        thrust_magnitude : float
//            Value of the constant thrust magnitude.
//        specific_impulse : float
//            Value of the constant specific impulse, used to link the thrust model to the mass propagation.
//        body_fixed_thrust_direction : numpy.ndarray, default=numpy.ndarray([])
//            Constant body-fixed thrust direction (positive x-direction by default). Note that this should be a unit-vector representing the direction opposite to the nozzle direction.
//
//        Returns
//        -------
//        ConstantThrustMagnitudeSettings
//            Constant thrust magnitude settings object.
//
//    )";
//
//
//    } else if(name == "custom_thrust_magnitude" && variant==1) {
//            return R"(
//
//        Create thrust magnitude settings from a custom thrust magnitude function.
//
//        Factory function that creates thrust magnitude from a custom thrust magnitude function.
//        This model defines a thrust force and specific impulse that can vary with time. The specific impulse is also
//        provided to apply a mass rate model in the propagation the vehicle dynamics, relating the thrust to the mass
//        decrease of the vehicle. Note that, if you wish to use a constant value for any or all of the first three
//        arguments, lambda expression can be used. Presently, the definition of the thrust direction in the body-fixed
//        frame is also defined through these derived classes. In essence, the ThrustMagnitudeSettings defines all local
//        (to the vehicle systems) settings for the thrust, while ThrustDirectionGuidanceSettings defines how the full
//        vehicle must orient itself in space for the required thrust direction to be achieved. At present, there is no
//        direct option for thrust-vector control (i.e. modifying the thrust direction in the body-fixed frame).
//
//
//        Parameters
//        ----------
//        thrust_magnitude_function : Callable[[float], float]
//            Function of time returning the value of the thrust magnitude.
//        specific_impulse_function : Callable[[float], float]
//            Function of time returning the value of the specific impulse, useful to link the mass propagation to the thrust model.
//        is_engine_on_function : Callable[[float], bool], default=lambda t: true
//            Function of time returning a boolean, denoting  whether the thrust should be engaged at all (e.g. thrust is 0 N if it returns false). It is useful to link the mass propagation to the thrust model.
//        body_fixed_thrust_direction
//            None
//        Callable[[], numpy.ndarray], default=lambda t: numpy.ndarray([])
//            Constant body-fixed thrust direction (positive x-direction by default). Note that this function should be a unit-vector representing the direction opposite to the nozzle direction. This setting can be used to incorporate thrust-vector control (TVC) into the thrust.
//        custom_thrust_reset_function : Callable[[float], ], default=lambda t: None
//            Function of time that updates any relevant aspects of the environment/system models, called before retrieving the thrust magnitude, specific impulse, and body-fixed thrust direction.
//
//        Returns
//        -------
//        FromFunctionThrustMagnitudeSettings
//            From function thrust magnitude settings object.
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//
//
//
//
//namespace dependent_variable {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//    } else if(name == "PropagationDependentVariables") {
//         return R"(
//
//        Enumeration of available propagation dependent variables.
//
//        Enumeration of propagation dependent variables supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "VariableSettings") {
//         return R"(
//
//        Functional base class to define settings for variables.
//
//     )";
//
//
//
//    } else if(name == "SingleDependentVariableSaveSettings") {
//         return R"(
//
//        `VariableSettings`-derived class to define settings for dependent variables that are to be saved during propagation.
//
//     )";
//
//
//
//    } else if(name == "SingleAccelerationDependentVariableSaveSettings") {
//         return R"(
//
//        `SingleDependentVariableSaveSettings`-derived class to save a single acceleration (norm or vector) during propagation.
//
//     )";
//
//
//
//
//    } else if(name == "create" && variant==0) {
//            return R"(
//
//        Function to create settings for a generic dependent variable.
//
//        Function to create settings for a dependent variable. It creates objects that calculate dependent variables
//        from the objects that define their settings. It is usually not relevant nor useful for the user.
//
//
//        Parameters
//        ----------
//        dependent_variable_list : list[SingleDependentVariableSaveSettings]
//            List of dependent variables to be saved.
//        print_variable_indices : bool, default=True
//            Whether the types of dependent variables to be saved should be printed on the terminal.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "mach_number" && variant==0) {
//            return R"(
//
//        Function to add the Mach number to the dependent variables to save.
//
//        Function to add the Mach number to the dependent variables to save. Requires an aerodynamic acceleration to be acting on the vehicle.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with atmosphere with respect to which the Mach number is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "altitude" && variant==0) {
//            return R"(
//
//        Function to add the altitude to the dependent variables to save.
//
//        Function to add the altitude to the dependent variables to save. It requires an aerodynamic acceleration to be acting on the vehicle and it depends on the central body's shape.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with atmosphere with respect to which the altitude is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "airspeed" && variant==0) {
//            return R"(
//
//        Function to add the airspeed to the dependent variables to save.
//
//        Function to add the airspeed to the dependent variables to save. Requires an aerodynamic acceleration to be acting on the vehicle.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with atmosphere with respect to which the airspeed is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "density" && variant==0) {
//            return R"(
//
//        Function to add the local density to the dependent variables to save.
//
//        Function to add the density (at position of body undergoing acceleration) to the dependent variables to save. Requires an aerodynamic acceleration to be acting on the vehicle.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        body_with_atmosphere : str
//            Body with atmosphere with respect to which the density is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "relative_speed" && variant==0) {
//            return R"(
//
//        Function to add the relative speed (norm of the velocity vector) to the dependent variables to save.
//
//        Function to add the relative speed (norm of the velocity vector) with respect to a second body to the dependent variables to save. The relative speed is computed between the bodies' centers of mass.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        relative_body : str
//            Body with respect to which the relative speed is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//    } else if(name == "relative_speed" && variant==1) {
//            return R"(
//
//        Function to add the relative position vector to the dependent variables to save.
//
//        Function to add the relative position vector with respect to a second body to the dependent variables to save. The relative position is computed between the bodies' centers of mass.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        relative_body : str
//            Body with respect to which the relative position is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "relative_speed" && variant==0) {
//            return R"(
//
//        Function to add the relative speed (norm of the velocity vector) to the dependent variables to save.
//
//        Function to add the relative speed (norm of the velocity vector) with respect to a second body to the dependent variables to save. The relative speed is computed between the bodies' centers of mass.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        relative_body : str
//            Body with respect to which the relative speed is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//    } else if(name == "relative_speed" && variant==1) {
//            return R"(
//
//        Function to add the relative position vector to the dependent variables to save.
//
//        Function to add the relative position vector with respect to a second body to the dependent variables to save. The relative position is computed between the bodies' centers of mass.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        relative_body : str
//            Body with respect to which the relative position is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "relative_distance" && variant==0) {
//            return R"(
//
//        Function to add the relative distance (norm of the position vector) to the dependent variables to save.
//
//        Function to add the relative distance (norm of the position vector) with respect to a second body to the dependent variables to save. The relative distance is computed between the bodies' centers of mass.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        relative_body : str
//            Body with respect to which the relative distance is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "relative_velocity" && variant==0) {
//            return R"(
//
//        Function to add the relative velocity vector to the dependent variables to save.
//
//        Function to add the relative velocity vector with respect to a second body to the dependent variables to save. The relative distance is computed between the bodies' centers of mass.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        relative_body : str
//            Body with respect to which the relative velocity is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "keplerian_state" && variant==0) {
//            return R"(
//
//        Function to add the Keplerian state to the dependent variables to save.
//
//        Function to add the Keplerian state to the dependent variables to save. The Keplerian state is returned in this order: 1: Semi-major Axis. 2: Eccentricity. 3: Inclination. 4: Argument of Periapsis. 5. Right Ascension of the Ascending Node. 6: True Anomaly.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the Keplerian state is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "modified_equinoctial_state" && variant==0) {
//            return R"(
//
//        Function to add the modified equinoctial state to the dependent variables to save.
//
//        Function to add the modified equinoctial state to the dependent variables to save. The value of the parameter I is automatically chosen as +1 or -1, depending on whether the inclination is smaller or larger than 90 degrees.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the modified equinoctial state is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "single_acceleration" && variant==0) {
//            return R"(
//
//        Function to add an acceleration vector to the dependent variables to save.
//
//        Function to add an acceleration vector to the dependent variables to save.
//
//        Parameters
//        ----------
//        acceleration_type : AvailableAcceleration
//            Acceleration type to be saved.
//        body_undergoing_acceleration : str
//            Body undergoing acceleration.
//        body_exerting_acceleration : str
//            Body exerting acceleration.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "single_acceleration_norm" && variant==0) {
//            return R"(
//
//        Function to add a scalar acceleration (norm of the acceleration vector) to the dependent variables to save.
//
//        Function to add a scalar acceleration (norm of the acceleration vector) to the dependent variables to save.
//
//        Parameters
//        ----------
//        acceleration_type : AvailableAcceleration
//            Acceleration type to be saved (see `AvailableAcceleration` enum).
//        body_undergoing_acceleration : str
//            Body undergoing acceleration.
//        body_exerting_acceleration : str
//            Body exerting acceleration.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "total_acceleration_norm" && variant==0) {
//            return R"(
//
//        Function to add the total scalar acceleration (norm of the vector) acting on a body to the dependent variables to save.
//
//        Function to add the total scalar acceleration (norm of the vector) acting on a body to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body undergoing acceleration.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "total_acceleration" && variant==0) {
//            return R"(
//
//        Function to add the total acceleration vector acting on a body to the dependent variables to save.
//
//        Function to add the total acceleration vector acting on a body to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body undergoing acceleration.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "single_torque_norm" && variant==0) {
//            return R"(
//
//        Function to add a single torque (norm of the torque vector) to the dependent variables to save.
//
//        Function to add a single torque (norm of the torque vector) to the dependent variables to save. The altitude depends on the shape of the central body.
//
//        Parameters
//        ----------
//        torque_type : AvailableTorque
//            Torque type to be saved.
//        body_undergoing_torque : str
//            Body undergoing torque.
//        body_exerting_torque : str
//            Body exerting torque.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "single_torque" && variant==0) {
//            return R"(
//
//        Function to add a single torque vector to the dependent variables to save.
//
//        Function to add a single torque vector to the dependent variables to save. The altitude depends on the shape of the central body.
//
//        Parameters
//        ----------
//        torque_type : AvailableTorque
//            Torque type to be saved.
//        body_undergoing_torque : str
//            Body undergoing torque.
//        body_exerting_torque : str
//            Body exerting torque.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "total_torque_norm" && variant==0) {
//            return R"(
//
//        Function to add the total torque (norm of the torque vector) to the dependent variables to save.
//
//        Function to add the total torque (norm of the torque vector) to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "total_torque" && variant==0) {
//            return R"(
//
//        Function to add the total torque vector to the dependent variables to save.
//
//        Function to add the total torque vector to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "spherical_harmonic_terms_acceleration" && variant==0) {
//            return R"(
//
//        Function to add a single term of the spherical harmonic acceleration vector to the dependent variables to save.
//
//        Function to add single term of the spherical harmonic acceleration vector to the dependent variables to save.
//
//        Parameters
//        ----------
//        body_undergoing_acceleration : str
//            Body undergoing acceleration.
//        body_exerting_acceleration : str
//            Body exerting acceleration.
//        component_indices : list[tuple]
//            Tuples of (degree, order) indicating the terms to save.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "spherical_harmonic_terms_acceleration_norm" && variant==0) {
//            return R"(
//
//        Function to add a single term of the spherical harmonic acceleration (norm of the vector) to the dependent variables to save.
//
//        Function to add single term of the spherical harmonic acceleration (norm of the vector) to the dependent variables to save.
//
//        Parameters
//        ----------
//        body_undergoing_acceleration : str
//            Body undergoing acceleration.
//        body_exerting_acceleration : str
//            Body exerting acceleration.
//        component_indices : list[tuple]
//            Tuples of (degree, order) indicating the terms to save.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "aerodynamic_force_coefficients" && variant==0) {
//            return R"(
//
//        Function to add the aerodynamic force coefficients to the dependent variables to save.
//
//        Function to add the aerodynamic force coefficients to the dependent variables to save. It requires an aerodynamic acceleration acting on the vehicle. The coefficients are returned in the following order: drag force, side force, lift force.
//
//        Parameters
//        ----------
//        body : str
//            Body undergoing acceleration.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "aerodynamic_moment_coefficients" && variant==0) {
//            return R"(
//
//        Function to add the aerodynamic moment coefficients to the dependent variables to save.
//
//        Function to add the aerodynamic force coefficients to the dependent variables to save. It requires an aerodynamic torque acting on the vehicle. The coefficients are returned in the following order: C_l, C_m, C_n (respectively about the X, Y, Z axes of the body-fixed frame, see Mooij 1994).
//
//        Parameters
//        ----------
//        body : str
//            Body undergoing acceleration.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "latitude" && variant==0) {
//            return R"(
//
//        Function to add the latitude to the dependent variables to save.
//
//        Function to add the latitude to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the latitude is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "geodetic_latitude" && variant==0) {
//            return R"(
//
//        Function to add the geodetic latitude to the dependent variables to save.
//
//        Function to add the geodetic latitude to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the geodetic latitude is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "longitude" && variant==0) {
//            return R"(
//
//        Function to add the longitude to the dependent variables to save.
//
//        Function to add the longitude to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the longitude is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "heading_angle" && variant==0) {
//            return R"(
//
//        Function to add the heading angle to the dependent variables to save.
//
//        Function to add the heading angle to the dependent variables to save. The heading angle is the angle between the X-axis of the vertical frame and the XZ-plane in the groundspeed-based trajectory frame (see Mooij, 1994).
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the heading angle is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "flight_path_angle" && variant==0) {
//            return R"(
//
//        Function to add the flight path angle to the dependent variables to save.
//
//        Function to add the flight path angle to the dependent variables to save. The flight path angle is the angle between the X-axis of the groundspeed-based trajectory frame and the local horizontal plane defined in the vertical reference frame (see Mooij, 1994).
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the flight path angle is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "angle_of_attack" && variant==0) {
//            return R"(
//
//        Function to add the angle of attack to the dependent variables to save.
//
//        Function to add the angle of attack angle to the dependent variables to save. The angle of attack is the angle between the X-axis of the body-fixed reference frame and the XY plane in the groundspeed-based aerodynamic frame (see Mooij, 1994).
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the angle of attack is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "sideslip_angle" && variant==0) {
//            return R"(
//
//        Function to add the sideslip angle to the dependent variables to save.
//
//        Function to add the sideslip angle to the dependent variables to save. The sideslip angle is ???
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the sideslip angle is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "bank_angle" && variant==0) {
//            return R"(
//
//        Function to add the bank angle to the dependent variables to save.
//
//        Function to add the bank angle to the dependent variables to save. The bank angle is ???
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the bank angle is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "radiation_pressure" && variant==0) {
//            return R"(
//
//        Function to add the radiation pressure to the dependent variables to save.
//
//        Function to add the radiation pressure to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        radiating_body : str
//            Radiating body.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "total_gravity_field_variation_acceleration" && variant==0) {
//            return R"(
//
//        Function to add the total gravity field variation acceleration to the dependent variables to save.
//
//        Function to add the total gravity field variation acceleration to the dependent variables to save. This function does not distinguish between different sources of variations of the gravity field. To select only one contribution, look for the single gravity field variation acceleration.
//
//        Parameters
//        ----------
//        body_undergoing_acceleration : str
//            Body whose dependent variable should be saved.
//        body_exerting_acceleration : str
//            Body exerting the acceleration.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "single_gravity_field_variation_acceleration" && variant==0) {
//            return R"(
//
//        Function to add a single gravity field variation acceleration to the dependent variables to save.
//
//        Function to add a single gravity field variation acceleration to the dependent variables to save. This function does distinguish between different sources of variations of the gravity field, but not between terms of the spherical harmonic expansion. To select specific combinations of order and degree, look for the single per term gravity field variation acceleration.
//
//        Parameters
//        ----------
//        body_undergoing_acceleration : str
//            Body whose dependent variable should be saved.
//        body_exerting_acceleration : str
//            Body exerting the acceleration.
//        deformation_type : str
//            Deformation type (see BodyDeformationTypes).
//        identifier : str, default=""
//            Identifier for the deformation type.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "single_per_terms_gravity_field_variation_acceleration" && variant==0) {
//            return R"(
//
//        Function to add specific spherical harmonic terms of a single gravity field variation acceleration to the dependent variables to save.
//
//        Function to add specific spherical harmonic terms of a single gravity field variation acceleration to the dependent variables to save.
//
//        Parameters
//        ----------
//        body_undergoing_acceleration : str
//            Body whose dependent variable should be saved.
//        body_exerting_acceleration : str
//            Body exerting the acceleration.
//        component_indices : list[tuple]
//            Tuples of (degree, order) indicating the terms to save.
//        deformation_type : str
//            Deformation type (see BodyDeformationTypes).
//        identifier : str, default=""
//            Identifier for the deformation type.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "body_fixed_airspeed_velocity" && variant==0) {
//            return R"(
//
//        Function to add the airspeed velocity vector to the dependent variables to save.
//
//        Function to add the airspeed velocity vector to the dependent variables to save. The airspeed velocity is expressed with respect to a central body and returned in a frame fixed to the same central body. It requires the central body to have an atmosphere.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the airspeed is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "body_fixed_groundspeed_velocity" && variant==0) {
//            return R"(
//
//        Function to add the groundspeed velocity vector to the dependent variables to save.
//
//        Function to add the groundspeed velocity vector to the dependent variables to save. The groundspeed velocity is expressed with respect to a central body and returned in a frame fixed to the same central body.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the groundspeed is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "inertial_to_body_fixed_rotation_frame" && variant==0) {
//            return R"(
//
//        Function to add the rotation matrix from the inertial RF to the body-fixed RF to the dependent variables to save.
//
//        Function to add the rotation matrix from the inertial RF to the body-fixed RF to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body of interest.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "lvlh_to_inertial_rotation_matrix" && variant==0) {
//            return R"(
//
//        Function to add the rotation matrix from the Local Vertical, Local Horizontal (LVLH) RF to the inertial RF to the dependent variables to save.
//
//        Function to add the rotation matrix from the Local Vertical, Local Horizontal (LVLH) RF to the inertial RF to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the groundspeed is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "inertial_to_body_fixed_313_euler_angles" && variant==0) {
//            return R"(
//
//        Function to add the rotation matrix from the inertial RF to the body-fixed RF to the dependent variables to save.
//
//        Function to add the rotation matrix from the inertial RF to the body-fixed RF to the dependent variables to save. It uses a 313-Euler angles representation.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "intermediate_aerodynamic_rotation_matrix_variable" && variant==0) {
//            return R"(
//
//        Function to add the rotation matrix from the a base aerodynamic RF to a target aerodynamic RF to the dependent variables to save.
//
//        Function to add the rotation matrix from the a base aerodynamic RF to a target aerodynamic RF to the dependent variables to save. The aerodynamic RFs are collected in the AerodynamicsReferenceFrames enum.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        base_frame : str
//            Base reference frame for the transformation.
//        target_frame : str
//            Target reference frame for the transformation.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "periapsis_altitude" && variant==0) {
//            return R"(
//
//        Function to add the altitude of periapsis to the dependent variables to save.
//
//        Function to add the altitude of periapsis to the dependent variables to save. The altitude depends on the shape of the central body.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the altitude of periapsis is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "control_surface_deflection" && variant==0) {
//            return R"(
//
//        Function to add the altitude of periapsis to the dependent variables to save.
//
//        Function to add the altitude of periapsis to the dependent variables to save. The altitude depends on the shape of the central body.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        control_surface : str
//            Control surface whose deflection should be saved.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "central_body_fixed_spherical_position" && variant==0) {
//            return R"(
//
//        Function to add the spherical, body-fixed position to the dependent variables to save.
//
//        Function to add the spherical position to the dependent variables to save. The spherical position is expressed in the central body's body-fixed RF.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the spherical, body-fixed is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "central_body_fixed_cartesian_position" && variant==0) {
//            return R"(
//
//        Function to add the cartesian, body-fixed position to the dependent variables to save.
//
//        Function to add the cartesian position to the dependent variables to save. The cartesian position is expressed in the central body's body-fixed RF.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the cartesian, body-fixed is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "body_mass" && variant==0) {
//            return R"(
//
//        Function to add the body mass to the dependent variables to save.
//
//        Function to add the body mass to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "radiation_pressure_coefficient" && variant==0) {
//            return R"(
//
//        Function to add the radiation pressure coefficient to the dependent variables to save.
//
//        Function to add the radiation pressure coefficient to the dependent variables to save.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        emitting_body : str
//            Emitting body.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "local_temperature" && variant==0) {
//            return R"(
//
//        Function to add the local temperature to the dependent variables to save.
//
//        Function to add the local temperature to the dependent variables to save (at position of body undergoing acceleration). It requires an aerodynamic acceleration to be acting on the body.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "local_dynamic_pressure" && variant==0) {
//            return R"(
//
//        Function to add the local temperature to the dependent variables to save.
//
//        Function to add the local temperature to the dependent variables to save (at position of body undergoing acceleration). It requires an aerodynamic acceleration to be acting on the body.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "local_aerodynamic_heat_rate" && variant==0) {
//            return R"(
//
//        Function to add the local aerodynamic heat rate to the dependent variables to save.
//
//        Function to add the local aerodynamic heat rate felt by the vehicle based on the current velocity and atmospheric conditions to the dependent variables to save (at position of body undergoing acceleration). It requires an aerodynamic acceleration to be acting on the body.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "local_aerodynamic_g_load" && variant==0) {
//            return R"(
//
//        Function to add the total aerodynamic G-load to the dependent variables to save.
//
//        Function to add the total aerodynamic G-load induced by the aerodynamic acceleration to the dependent variables to save (at position of body undergoing acceleration). It requires an aerodynamic acceleration to be acting on the body.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "stagnation_point_heat_flux" && variant==0) {
//            return R"(
//
//        Function to add the heat flux at the stagnation point to the dependent variables to save.
//
//        Function to add the heat flux induced by atmospheric friction at the stagnation point to the dependent variables to save. It requires an aerodynamic acceleration to be acting on the body and a vehicle nose radius to be defined.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "total_mass_rate" && variant==0) {
//            return R"(
//
//        Function to add the total mass rate to the dependent variables to save.
//
//        Function to add the total mass rate to the dependent variables to save. It requires the body mass to be numerically propagated.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "aerodynamic_g_load" && variant==0) {
//            return R"(
//
//        Function to add the aerodynamic g-load to the dependent variables to save.
//
//        Function to add the aerodynamic g-load to the dependent variables to save (at position of body undergoing acceleration). It requires an aerodynamic acceleration to be acting on the body.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the aerodynamic g-load is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "dynamic_pressure" && variant==0) {
//            return R"(
//
//        Function to add the dynamic pressure to the dependent variables to save.
//
//        Function to add the dynamic pressure to the dependent variables to save. It requires the central body to have an atmosphere.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the dynamic pressure is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else if(name == "atmospheric_temperature" && variant==0) {
//            return R"(
//
//        Function to add the atmospheric temperature to the dependent variables to save.
//
//        Function to add the atmospheric temperature to the dependent variables to save. It requires the central body to have an atmosphere.
//
//        Parameters
//        ----------
//        body : str
//            Body whose dependent variable should be saved.
//        central_body : str
//            Body with respect to which the atmospheric temperature is computed.
//
//        Returns
//        -------
//        SingleDependentVariableSaveSettings
//            Dependent variable settings object.
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//
//
//
//
//namespace integrator {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//    } else if(name == "AvailableIntegrators") {
//         return R"(
//
//        Enumeration of available integrators.
//
//        Enumeration of integrators supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "CoefficientSets") {
//         return R"(
//
//        Coefficient sets for Runge-Kutta integrators.
//
//        Coefficient sets for Runge-Kutta integrators.
//
//     )";
//
//
//
//    } else if(name == "ExtrapolationMethodStepSequences") {
//         return R"(
//
//        Enumeration of available extrapolation method step sequences.
//
//        Enumeration of extrapolation method step sequences supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "IntegratorSettings") {
//         return R"(
//
//        Functional base class to define settings for integrators.
//
//     )";
//
//
//
//    } else if(name == "RungeKuttaVariableStepSizeSettingsScalarTolerances") {
//         return R"(
//
//        `IntegratorSettings`-derived class to define settings for Runge Kutta integrators with scalar tolerances.
//
//     )";
//
//
//
//    } else if(name == "RungeKuttaVariableStepSizeSettingsVectorTolerances") {
//         return R"(
//
//        `IntegratorSettings`-derived class to define settings for Runge Kutta integrators with vector tolerances.
//
//     )";
//
//
//
//    } else if(name == "BulirschStoerIntegratorSettings") {
//         return R"(
//
//        `IntegratorSettings`-derived class to define settings for Bulirsch-Stoer integrator settings.
//
//     )";
//
//
//
//    } else if(name == "AdamsBashforthMoultonSettings") {
//         return R"(
//
//        `IntegratorSettings`-derived class to define settings for Adams-Bashforth-Moulton integrator settings.
//
//     )";
//
//
//
//
//    } else if(name == "euler" && variant==0) {
//            return R"(
//
//        Creates the settings for the Euler integrator.
//
//        Factory function to create settings for the Euler integrator. For this integrator, the step size is kept
//        constant.
//
//
//        Parameters
//        ----------
//        initial_time : float
//            Start time (independent variable) of numerical integration.
//        initial_time_step : float
//            Initial and constant value for the time step.
//        save_frequency : int, default=1
//            Frequency at which to save the numerical integrated states (expressed per unit integration time step, with n = saveFrequency, so n = 1 means that the state is saved once per integration step).
//        assess_termination_on_minor_steps : bool, default=false
//            Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
//
//        Returns
//        -------
//        IntegratorSettings
//            Integrator settings object.
//
//    )";
//
//
//
//    } else if(name == "runge_kutta_4" && variant==0) {
//            return R"(
//
//        Creates the settings for the Runge Kutta 4 integrator.
//
//        Factory function to create settings for the Runge Kutta 4 integrator. For this integrator, the step size is kept
//        constant.
//
//
//        Parameters
//        ----------
//        initial_time : float
//            Start time (independent variable) of numerical integration.
//        initial_time_step : float
//            Initial and constant value for the time step.
//        save_frequency : int, default=1
//            Frequency at which to save the numerical integrated states (expressed per unit integration time step, with n = saveFrequency, so n = 1 means that the state is saved once per integration step).
//        assess_termination_on_minor_steps : bool, default=false
//            Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
//
//        Returns
//        -------
//        IntegratorSettings
//            Integrator settings object.
//
//    )";
//
//
//
//    } else if(name == "rungeKuttaVariableStepSettingsScalarTolerances" && variant==0) {
//            return R"(
//
//    )";
//
//
//
//    } else if(name == "rungeKuttaVariableStepSettingsVectorTolerances" && variant==0) {
//            return R"(
//
//    )";
//
//
//
//    } else if(name == "runge_kutta_variable_step_size" && variant==0) {
//            return R"(
//
//        Creates the settings for the Runge-Kutta variable step size integrator.
//
//        Factory function to create settings for the Runge-Kutta variable step size integrator with vector tolerances. # [py]
//        For this integrator, the step size is varied based on the tolerances and safety factor provided. # [py]
//        The tolerance can be either scalar or vector; it is composed of an absolute and a relative part. # [py]
//        Different coefficient sets (Butcher's tableau) can be used (see the `RungeKuttaCoefficients::CoefficientSets` enum). # [py]
//
//
//        Parameters
//        ----------
//        initial_time : float
//            Start time (independent variable) of numerical integration.
//        initial_time_step : float
//            Initial time step to be used.
//        coefficient_set : RungeKuttaCoefficients::CoefficientSets
//            Coefficient set (Butcher's tableau) to be used in the integration.
//        minimum_step_size : float
//            Minimum time step to be used during the integration.
//        maximum_step_size : float
//            Maximum time step to be used during the integration.
//        relative_error_tolerance : float or np.ndarray
//            Relative vector tolerance to adjust the time step.
//        absolute_error_tolerance : float or np.ndarray
//            Absolute vector tolerance to adjust the time step.
//        save_frequency : int, default=1
//            Frequency at which to save the numerical integrated states (expressed per unit integration time step, with n = saveFrequency, so n = 1 means that the state is saved once per integration step).
//        assess_termination_on_minor_steps : bool, default=false
//            Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
//        safety_factor : float, default=0.8
//            Safety factor used in the step size control.
//        maximum_factor_increase : float, default=4.0
//            Maximum increase between consecutive time steps, expressed as the factor between new and old step size.
//        minimum_factor_increase : float, default=0.1
//            Minimum increase between consecutive time steps, expressed as the factor between new and old step size.
//
//        Returns
//        -------
//        RungeKuttaVariableStepSizeSettingsScalarTolerances or RungeKuttaVariableStepSizeSettingsVectorTolerances
//            RungeKuttaVariableStepSizeSettingsScalarTolerances or RungeKuttaVariableStepSizeSettingsVectorTolerances object.
//
//    )";
//
//
//
//    } else if(name == "bulirsch_stoer" && variant==0) {
//            return R"(
//
//        Creates the settings for the Bulirsch-Stoer integrator.
//
//        Factory function to create settings for the Bulirsch-Stoer integrator.
//        For this integrator, the step size is varied based on the tolerances and safety factor provided.
//        The tolerance is composed of an absolute and a relative part.
//        Different extrapolation sequences can be used (see the `ExtrapolationMethodStepSequences` enum).
//
//
//        Parameters
//        ----------
//        initial_time : float
//            Start time (independent variable) of numerical integration.
//        initial_time_step : float
//            Initial time step to be used.
//        extrapolation_sequence : ExtrapolationMethodStepSequences
//            Extrapolation sequence to be used in the integration.
//        maximum_number_of_steps : int
//            Number of entries in the sequence (e.g., number of integrations used for a single extrapolation).
//        minimum_step_size : float
//            Minimum time step to be used during the integration.
//        maximum_step_size : float
//            Maximum time step to be used during the integration.
//        relative_error_tolerance : float, default=1.0E-12
//            Relative tolerance to adjust the time step.
//        absolute_error_tolerance : float, default=1.0E-12
//            Relative tolerance to adjust the time step.
//        save_frequency : int, default=1
//            Frequency at which to save the numerical integrated states (expressed per unit integration time step, with n = saveFrequency, so n = 1 means that the state is saved once per integration step).
//        assess_termination_on_minor_steps : bool, default=false
//            Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
//        safety_factor : float, default=0.7
//            Safety factor used in the step size control.
//        maximum_factor_increase : float, default=10.0
//            Maximum increase between consecutive time steps, expressed as the factor between new and old step size.
//        minimum_factor_increase : float, default=0.1
//            Minimum increase between consecutive time steps, expressed as the factor between new and old step size.
//
//        Returns
//        -------
//        BulirschStoerIntegratorSettings
//            BulirschStoerIntegratorSettings object.
//
//    )";
//
//
//
//    } else if(name == "adams_bashforth_moulton" && variant==0) {
//            return R"(
//
//        Creates the settings for the Bulirsch-Stoer integrator.
//
//        Factory function to create settings for the Adams-Bashorth-Moulton integrator.
//        For this integrator, the step size is varied based on the tolerances and safety factor provided.
//        The tolerance is composed of an absolute and a relative part.
//        Different coefficient sets (Butcher's tableau) can be used (see the `RungeKuttaCoefficients::CoefficientSets` enum).
//
//
//        Parameters
//        ----------
//        initial_time : float
//            Start time (independent variable) of numerical integration.
//        initial_time_step : float
//            Initial time step to be used.
//        minimum_step_size : float
//            Minimum time step to be used during the integration.
//        maximum_step_size : float
//            Maximum time step to be used during the integration.
//        relative_error_tolerance : float, default=1.0E-12
//            Relative tolerance to adjust the time step.
//        absolute_error_tolerance : float, default=1.0E-12
//            Relative tolerance to adjust the time step.
//        minimum_order
//            Minimum order of the integrator.
//        maximum_order
//            Maximum order of the integrator.
//        save_frequency : int, default=1
//            Frequency at which to save the numerical integrated states (expressed per unit integration time step, with n = saveFrequency, so n = 1 means that the state is saved once per integration step).
//        assess_termination_on_minor_steps : bool, default=false
//            Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
//        bandwidth : float, default=200.0
//            Maximum error factor for doubling the stepsize.
//
//        Returns
//        -------
//        AdamsBashforthMoultonSettings
//            AdamsBashforthMoultonSettings object.
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//
//
//
//
//namespace mass_rate {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//    } else if(name == "AvailableMassRateModels") {
//         return R"(
//
//        Enumeration of available mass rate models.
//
//        Enumeration of mass rate models supported by tuday.
//
//     )";
//
//
//
//    } else if(name == "MassRateModelSettings") {
//         return R"(
//
//        Functional base class to define settings for mass rates.
//
//     )";
//
//
//
//    } else if(name == "FromThrustMassRateSettings") {
//         return R"(
//
//        `MassRateModelSettings`-derived class to define settings for a mass rate model derived from a thrust model.
//
//     )";
//
//
//
//    } else if(name == "CustomMassRateSettings") {
//         return R"(
//
//        `MassRateModelSettings`-derived class to define settings for a custom mass rate model.
//
//     )";
//
//
//
//
//    } else if(name == "from_thrust" && variant==0) {
//            return R"(
//
//        Creates the settings for a mass rate model defined from a thrust model.
//
//        Creates the settings for a mass rate model defined from a thrust model. The mass rate model is derived from
//        the associated body's engine model. It is possible to consider only a specific engine or all engines.
//
//
//        Parameters
//        ----------
//        use_all_thrust_models : bool, default=true
//            Denotes whether all engines of the associated body are to be combined into a single thrust model.
//        associated_thrust_source : str, default=""
//            Name of engine model from which thrust is to be derived (must be empty if the first argument is set to true).
//
//        Returns
//        -------
//        FromThrustMassRateSettings
//            From thrust mass rate settings object.
//
//    )";
//
//
//
//    } else if(name == "custom" && variant==0) {
//            return R"(
//
//        Creates the settings for a mass rate model defined from a thrust model.
//
//        Creates the settings for a custom mass rate model defined through a mass rate function. The function must have
//        time as an independent variable.
//
//
//        Parameters
//        ----------
//        mass_rate_function : Callable[[float], float]
//            Function of time defining the custom mass rate.
//
//        Returns
//        -------
//        CustomMassRateSettings
//            Custom mass rate settings object.
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//
//
//
//
//namespace propagator {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//    } else if(name == "TranslationalPropagatorType") {
//         return R"(
//
//        Enumeration of available translational propagator types.
//
//        Enumeration of translational propagator types supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "RotationalPropagatorType") {
//         return R"(
//
//        Enumeration of available rotational propagator types.
//
//        Enumeration of rotational propagator types supported by tudat. # [py]
//
//     )";
//
//
//
//    } else if(name == "StateType") {
//         return R"(
//
//        Enumeration of available integrated state types.
//
//        Enumeration of integrated state types supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "RotationalPropagatorType") {
//         return R"(
//
//        Enumeration of available rotational propagator types.
//
//        Enumeration of rotational propagator types supported by tudat. # [py]
//
//     )";
//
//
//
//    } else if(name == "DependentVariableSaveSettings") {
//         return R"(
//
//        Functional class to define settings for dependent variable to save.
//
//     )";
//
//
//    } else if(name == "DependentVariableSaveSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user because this is a base class.
//
//    )";
//
//
//
//
//    } else if(name == "PropagatorSettings") {
//         return R"(
//
//        Functional base class to define settings for propagators.
//
//     )";
//
//
//    } else if(name == "PropagatorSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user because this is a base class.
//
//    )";
//
//
//
//    } else if(name == "PropagatorSettings.reset_initial_states" && variant==0) {
//            return R"(
//
//        Function to reset the initial state used as input for numerical integration.
//
//        Function to reset the initial state used as input for numerical integration.
//
//
//        Parameters
//        ----------
//        initial_states : numpy.ndarray
//            Initial states to be reset for the numerical propagation.
//    )";
//
//
//
//
//    } else if(name == "MultiArcPropagatorSettings") {
//         return R"(
//
//        `PropagatorSettings`-derived class to define settings for multi-arc dynamics.
//
//     )";
//
//
//    } else if(name == "MultiArcPropagatorSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user. Settings objects for integrators should be
//        instantiated through the factory functions of a derived class.
//
//    )";
//
//
//
//
//    } else if(name == "HybridArcPropagatorSettings") {
//         return R"(
//
//        `PropagatorSettings`-derived class to define settings for hybrid-arc dynamics.
//
//     )";
//
//
//    } else if(name == "HybridArcPropagatorSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user. Settings objects for integrators should be
//        instantiated through the factory functions of a derived class.
//
//    )";
//
//
//
//
//    } else if(name == "SingleArcPropagatorSettings") {
//         return R"(
//
//        `PropagatorSettings`-derived class to define settings for single-arc dynamics.
//
//     )";
//
//
//    } else if(name == "SingleArcPropagatorSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user. Settings objects for integrators should be
//        instantiated through the factory functions of a derived class.
//
//    )";
//
//
//
//
//    } else if(name == "TranslationalStatePropagatorSettings") {
//         return R"(
//
//        `SingleArcPropagatorSettings`-derived class to define settings for single-arc translational dynamics.
//
//     )";
//
//
//    } else if(name == "TranslationalStatePropagatorSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user. Settings objects for integrators should be
//        instantiated through the factory functions of a derived class.
//
//    )";
//
//
//
//    } else if(name == "TranslationalStatePropagatorSettings.reset_initial_states" && variant==0) {
//            return R"(
//
//        Function to reset the initial state used as input for numerical integration.
//
//        Function to reset the initial state used as input for numerical integration.
//
//
//        Parameters
//        ----------
//        initial_states : numpy.ndarray
//            Initial states to be reset for the numerical propagation.
//    )";
//
//
//
//    } else if(name == "TranslationalStatePropagatorSettings.recreate_state_derivative_models" && variant==0) {
//            return R"(
//
//        Function to (re)create the integrated state models (e.g. acceleration/torque/mass models).
//
//        Function to create the integrated state models (e.g. acceleration/torque/mass models) for
//        each fo the propagators state types contained in `propagatorSettingsMap_`.
//
//
//        Parameters
//        ----------
//        bodies : SystemOfBodies
//            System of bodies used in the propagation.
//    )";
//
//
//
//    } else if(name == "TranslationalStatePropagatorSettings.single_type_settings" && variant==0) {
//            return R"(
//
//        Function to retrieve a single type of propagator.
//
//        Function to retrieve a single type of propagator (translational, rotational or mass). This function is
//        often used in multi-type propagation.
//
//
//        Parameters
//        ----------
//        state_type : IntegratedStateType
//            State type to be retrieved.
//    )";
//
//
//
//
//    } else if(name == "RotationalStatePropagatorSettings") {
//         return R"(
//
//        `SingleArcPropagatorSettings`-derived class to define settings for single-arc rotational state propagation.
//
//     )";
//
//
//    } else if(name == "RotationalStatePropagatorSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user. Settings objects for integrators should be
//        instantiated through the factory functions of a derived class.
//
//    )";
//
//
//
//
//    } else if(name == "MultiTypePropagatorSettings") {
//         return R"(
//
//        `SingleArcPropagatorSettings`-derived class to define settings for propagation of multiple quantities.
//
//     )";
//
//
//    } else if(name == "MultiTypePropagatorSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user. Settings objects for integrators should be
//        instantiated through the factory functions of a derived class.
//
//    )";
//
//
//
//    } else if(name == "MultiTypePropagatorSettings.reset_initial_states" && variant==0) {
//            return R"(
//
//        Function to reset the initial states used as input for numerical integration.
//
//        Function to reset the initial states used as input for numerical integration.
//
//
//        Parameters
//        ----------
//        initial_states : numpy.ndarray
//            Initial states to be reset for the numerical propagation.
//    )";
//
//
//
//    } else if(name == "MultiTypePropagatorSettings.recreate_state_derivative_models" && variant==0) {
//            return R"(
//
//        Function to (re)create the integrated state models (e.g. acceleration/torque/mass models).
//
//        Function to create the integrated state models (e.g. acceleration/torque/mass models) for
//        each of the propagators state types contained in `propagatorSettingsMap_`.
//
//
//        Parameters
//        ----------
//        bodies : SystemOfBodies
//            System of bodies used in the propagation.
//    )";
//
//
//
//    } else if(name == "MultiTypePropagatorSettings.single_type_settings" && variant==0) {
//            return R"(
//
//        Function to retrieve a single type of propagator.
//
//        Function to retrieve a single type of propagator (translational, rotational or mass). This function is
//        often used in multi-type propagation.
//
//
//        Parameters
//        ----------
//        state_type : IntegratedStateType
//            State type to be retrieved.
//    )";
//
//
//
//
//    } else if(name == "PropagationTerminationSettings") {
//         return R"(
//
//        Functional base class to define termination settings for the propagation.
//
//     )";
//
//
//    } else if(name == "PropagationTerminationSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user because this is a base class.
//
//    )";
//
//
//
//
//    } else if(name == "PropagationDependentVariableTerminationSettings") {
//         return R"(
//
//        `PropagationTerminationSettings`-derived class to define termination settings for the propagation from dependent variables.
//
//     )";
//
//
//    } else if(name == "PropagationDependentVariableTerminationSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user. Settings objects for integrators should be
//        instantiated through the factory functions of a derived class.
//
//    )";
//
//
//
//
//    } else if(name == "PropagationTimeTerminationSettings") {
//         return R"(
//
//        `PropagationTerminationSettings`-derived class to define termination settings for the propagation from propagation time.
//
//     )";
//
//
//    } else if(name == "PropagationTimeTerminationSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user. Settings objects for integrators should be
//        instantiated through the factory functions of a derived class.
//
//    )";
//
//
//
//
//    } else if(name == "PropagationCPUTimeTerminationSettings") {
//         return R"(
//
//        `PropagationTerminationSettings`-derived class to define termination settings for the propagation from CPU time.
//
//     )";
//
//
//    } else if(name == "PropagationCPUTimeTerminationSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user. Settings objects for integrators should be
//        instantiated through the factory functions of a derived class.
//
//    )";
//
//
//
//
//    } else if(name == "PropagationCustomTerminationSettings") {
//         return R"(
//
//        `PropagationTerminationSettings`-derived class to define custom termination settings for the propagation.
//
//     )";
//
//
//    } else if(name == "PropagationCustomTerminationSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user. Settings objects for integrators should be
//        instantiated through the factory functions of a derived class.
//
//    )";
//
//
//
//
//    } else if(name == "PropagationHybridTerminationSettings") {
//         return R"(
//
//        `PropagationTerminationSettings`-derived class to define hybrid termination settings for the propagation.
//
//     )";
//
//
//    } else if(name == "PropagationHybridTerminationSettings.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user. Settings objects for integrators should be
//        instantiated through the factory functions of a derived class.
//
//    )";
//
//
//
//
//    } else if(name == "SingleArcDynamicsSimulator") {
//         return R"(
//
//        Base class to define single arc dynamics simulator settings.
//
//     )";
//
//
//    } else if(name == "SingleArcDynamicsSimulator.ctor" && variant==0) {
//            return R"(
//
//        Constructor.
//
//        Instances of this class are typically not generated by the user because this is a base class.
//
//    )";
//
//
//
//    } else if(name == "SingleArcDynamicsSimulator.integrate_equations_of_motion" && variant==0) {
//            return R"(
//
//        This function numerically (re-)integrates the equations of motion.
//
//        This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
//        and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
//        attribute.
//
//
//        Parameters
//        ----------
//        initial_states : numpy.ndarray
//            Initial state vector that is to be used for numerical integration. Note that this state should
//            be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_), but not in the propagator-
//            specific form (i.e Encke, Gauss, etc. for translational dynamics).
//
//    )";
//
//
//    } else if(name == "SingleArcDynamicsSimulator.integrate_equations_of_motion" && variant==1) {
//            return R"(
//
//        This function numerically (re-)integrates the equations of motion.
//
//        This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
//        and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
//        attribute.
//
//
//        Parameters
//        ----------
//        initial_states : numpy.ndarray
//            Initial state vector that is to be used for numerical integration. Note that this state should
//            be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_), but not in the propagator-
//            specific form (i.e Encke, Gauss, etc. for translational dynamics).
//
//    )";
//
//
//
//    } else if(name == "SingleArcDynamicsSimulator.integrate_equations_of_motion" && variant==0) {
//            return R"(
//
//        This function numerically (re-)integrates the equations of motion.
//
//        This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
//        and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
//        attribute.
//
//
//        Parameters
//        ----------
//        initial_states : numpy.ndarray
//            Initial state vector that is to be used for numerical integration. Note that this state should
//            be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_), but not in the propagator-
//            specific form (i.e Encke, Gauss, etc. for translational dynamics).
//
//    )";
//
//
//    } else if(name == "SingleArcDynamicsSimulator.integrate_equations_of_motion" && variant==1) {
//            return R"(
//
//        This function numerically (re-)integrates the equations of motion.
//
//        This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
//        and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
//        attribute.
//
//
//        Parameters
//        ----------
//        initial_states : numpy.ndarray
//            Initial state vector that is to be used for numerical integration. Note that this state should
//            be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_), but not in the propagator-
//            specific form (i.e Encke, Gauss, etc. for translational dynamics).
//
//    )";
//
//
//
//
//
//    } else if(name == "combine_initial_states" && variant==0) {
//            return R"(
//
//        Function to retrieve the initial state for a list of propagator settings.
//
//        Function to retrieve the initial state for a list of propagator settings. This way, the initial state for
//        different quantities to be propagated (e.g., translational state, rotational state, mass) are retrieved and
//        organized in a single container.
//
//
//        Parameters
//        ----------
//        propagator_settings_per_type : dict
//            Propagator settings where the type of propagation is reported as key and the respective list of propagator settings as value.
//
//        Returns
//        -------
//        numpy.ndarray
//            Vector of initial states, sorted in order of IntegratedStateType, and then in the order of the vector of SingleArcPropagatorSettings of given type.
//
//    )";
//
//
//
//    } else if(name == "create_acceleration_models" && variant==0) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies and central bodies are provided through a dictionary.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_acceleration_per_body : SelectedAccelerationMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        central_bodies : dict
//            Key-value container indicating the body to propagate as key and its central body as value.
//
//        Returns
//        -------
//        AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==1) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies and central bodies are provided as two separate lists with the same order.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_acceleration_per_body : SelectedAccelerationMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        bodies_to_propagate : list
//            List of bodies to propagate.
//        central_bodies : list
//            List of central bodies, each referred to each propagated body in the same order.
//
//        Returns
//        -------
//        AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==2) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies is provided as a list.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_torque_per_body : SelectedTorqueMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        bodies_to_propagate : list
//            List of bodies to propagate.
//
//        Returns
//        -------
//        TorqueModelMap
//            Set of torques acting on the bodies to propagate, provided as torque models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==3) {
//            return R"(
//
//        Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.
//
//        Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.
//
//
//        Parameters
//        ----------
//        bodies_to_propagate : List[str]
//            List of bodies to be propagated.
//        central_bodies : List[str]
//            List of central bodies, each referred to a body being propagated (in the same order).
//        bodies_to_propagate : body_system
//            System of bodies used in the propagation.
//        initial_time : float
//            Initial time of the propagation.
//
//        Returns
//        -------
//        numpy.ndarray
//            Time at which the states should be retrieved.
//
//    )";
//
//
//
//    } else if(name == "create_acceleration_models" && variant==0) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies and central bodies are provided through a dictionary.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_acceleration_per_body : SelectedAccelerationMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        central_bodies : dict
//            Key-value container indicating the body to propagate as key and its central body as value.
//
//        Returns
//        -------
//        AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==1) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies and central bodies are provided as two separate lists with the same order.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_acceleration_per_body : SelectedAccelerationMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        bodies_to_propagate : list
//            List of bodies to propagate.
//        central_bodies : list
//            List of central bodies, each referred to each propagated body in the same order.
//
//        Returns
//        -------
//        AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==2) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies is provided as a list.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_torque_per_body : SelectedTorqueMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        bodies_to_propagate : list
//            List of bodies to propagate.
//
//        Returns
//        -------
//        TorqueModelMap
//            Set of torques acting on the bodies to propagate, provided as torque models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==3) {
//            return R"(
//
//        Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.
//
//        Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.
//
//
//        Parameters
//        ----------
//        bodies_to_propagate : List[str]
//            List of bodies to be propagated.
//        central_bodies : List[str]
//            List of central bodies, each referred to a body being propagated (in the same order).
//        bodies_to_propagate : body_system
//            System of bodies used in the propagation.
//        initial_time : float
//            Initial time of the propagation.
//
//        Returns
//        -------
//        numpy.ndarray
//            Time at which the states should be retrieved.
//
//    )";
//
//
//
//    } else if(name == "create_acceleration_models" && variant==0) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies and central bodies are provided through a dictionary.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_acceleration_per_body : SelectedAccelerationMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        central_bodies : dict
//            Key-value container indicating the body to propagate as key and its central body as value.
//
//        Returns
//        -------
//        AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==1) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies and central bodies are provided as two separate lists with the same order.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_acceleration_per_body : SelectedAccelerationMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        bodies_to_propagate : list
//            List of bodies to propagate.
//        central_bodies : list
//            List of central bodies, each referred to each propagated body in the same order.
//
//        Returns
//        -------
//        AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==2) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies is provided as a list.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_torque_per_body : SelectedTorqueMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        bodies_to_propagate : list
//            List of bodies to propagate.
//
//        Returns
//        -------
//        TorqueModelMap
//            Set of torques acting on the bodies to propagate, provided as torque models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==3) {
//            return R"(
//
//        Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.
//
//        Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.
//
//
//        Parameters
//        ----------
//        bodies_to_propagate : List[str]
//            List of bodies to be propagated.
//        central_bodies : List[str]
//            List of central bodies, each referred to a body being propagated (in the same order).
//        bodies_to_propagate : body_system
//            System of bodies used in the propagation.
//        initial_time : float
//            Initial time of the propagation.
//
//        Returns
//        -------
//        numpy.ndarray
//            Time at which the states should be retrieved.
//
//    )";
//
//
//
//    } else if(name == "create_acceleration_models" && variant==0) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies and central bodies are provided through a dictionary.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_acceleration_per_body : SelectedAccelerationMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        central_bodies : dict
//            Key-value container indicating the body to propagate as key and its central body as value.
//
//        Returns
//        -------
//        AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==1) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies and central bodies are provided as two separate lists with the same order.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_acceleration_per_body : SelectedAccelerationMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        bodies_to_propagate : list
//            List of bodies to propagate.
//        central_bodies : list
//            List of central bodies, each referred to each propagated body in the same order.
//
//        Returns
//        -------
//        AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==2) {
//            return R"(
//
//        Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
//
//        Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
//        bodies is provided as a list.
//
//
//        Parameters
//        ----------
//        body_system : SystemOfBodies
//            System of bodies to be used in the propagation.
//        selected_torque_per_body : SelectedTorqueMap
//            Key-value container indicating the acceleration type(s) as value and the body undergoing such acceleration(s) as key.
//        bodies_to_propagate : list
//            List of bodies to propagate.
//
//        Returns
//        -------
//        TorqueModelMap
//            Set of torques acting on the bodies to propagate, provided as torque models.
//
//    )";
//
//
//    } else if(name == "create_acceleration_models" && variant==3) {
//            return R"(
//
//        Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.
//
//        Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.
//
//
//        Parameters
//        ----------
//        bodies_to_propagate : List[str]
//            List of bodies to be propagated.
//        central_bodies : List[str]
//            List of central bodies, each referred to a body being propagated (in the same order).
//        bodies_to_propagate : body_system
//            System of bodies used in the propagation.
//        initial_time : float
//            Initial time of the propagation.
//
//        Returns
//        -------
//        numpy.ndarray
//            Time at which the states should be retrieved.
//
//    )";
//
//
//
//    } else if(name == "translational" && variant==0) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==1) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved, provided as a list of SingleDependentVariableSaveSettings objects (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==2) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==3) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==4) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==5) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==6) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "translational" && variant==0) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==1) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved, provided as a list of SingleDependentVariableSaveSettings objects (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==2) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==3) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==4) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==5) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==6) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "translational" && variant==0) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==1) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved, provided as a list of SingleDependentVariableSaveSettings objects (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==2) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==3) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==4) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==5) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==6) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "translational" && variant==0) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==1) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved, provided as a list of SingleDependentVariableSaveSettings objects (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==2) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==3) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==4) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==5) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==6) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "translational" && variant==0) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==1) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved, provided as a list of SingleDependentVariableSaveSettings objects (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==2) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==3) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==4) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==5) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==6) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "translational" && variant==0) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==1) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved, provided as a list of SingleDependentVariableSaveSettings objects (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==2) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==3) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==4) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==5) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==6) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "translational" && variant==0) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==1) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved, provided as a list of SingleDependentVariableSaveSettings objects (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==2) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==3) {
//            return R"(
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==4) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==5) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "translational" && variant==6) {
//            return R"(
//
//        Factory function to create translational state propagator settings with time termination conditions.
//
//        Factory function to create translational state propagator settings with generic stopping conditions.
//        It works by providing a key-value acceleration container, containing the list of accelerations acting on
//        each body. The map has as key a string denoting the name of the body on which a set of accelerations, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the acceleration
//        and the value (a pointer to) an acceleration settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the final time provided.
//
//
//        Parameters
//        ----------
//        central_bodies : list[str]
//            List of central bodies with respect to which the bodies to be integrated are propagated.
//        acceleration_models : AccelerationMap
//            Set of accelerations acting on the bodies to propagate, provided as acceleration settings objects.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_time : float
//            Final time of the propagation to be used as termination criterion.
//        propagator : TranslationalPropagatorType, default=cowell
//            Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        TranslationalStatePropagatorSettings
//            Translational state propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "rotational" && variant==0) {
//            return R"(
//
//        Factory function to create rotational state propagator settings with generic stopping conditions.
//
//        Factory function to create rotational state propagator settings with generic stopping conditions.
//        It works by providing a key-value torque container, containing the list of torques acting on
//        each body. The map has as key a string denoting the name of the body on which a set of torques, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the torque
//        and the value (a pointer to) a torque model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        torque_models : TorqueModelMap
//            Set of torques acting on the bodies to propagate, provided as torque models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : RotationalPropagatorType, default=quaternions
//            Type of rotational propagator to be used (see `RotationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        RotationalStatePropagatorSettings
//            Rotational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "rotational" && variant==1) {
//            return R"(
//
//        Factory function to create rotational state propagator settings with generic stopping conditions.
//
//        Factory function to create rotational state propagator settings with generic stopping conditions.
//        It works by providing a key-value torque container, containing the list of torques acting on
//        each body. The map has as key a string denoting the name of the body on which a set of torques, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the torque
//        and the value (a pointer to) a torque settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        torque_settings : SelectedTorqueMap
//            Set of torques acting on the bodies to propagate, provided as torque settings object.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : RotationalPropagatorType, default=quaternions
//            Type of rotational propagator to be used (see `RotationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        RotationalStatePropagatorSettings
//            Rotational state propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "rotational" && variant==0) {
//            return R"(
//
//        Factory function to create rotational state propagator settings with generic stopping conditions.
//
//        Factory function to create rotational state propagator settings with generic stopping conditions.
//        It works by providing a key-value torque container, containing the list of torques acting on
//        each body. The map has as key a string denoting the name of the body on which a set of torques, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the torque
//        and the value (a pointer to) a torque model. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        torque_models : TorqueModelMap
//            Set of torques acting on the bodies to propagate, provided as torque models.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : RotationalPropagatorType, default=quaternions
//            Type of rotational propagator to be used (see `RotationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        RotationalStatePropagatorSettings
//            Rotational state propagator settings object.
//
//    )";
//
//
//    } else if(name == "rotational" && variant==1) {
//            return R"(
//
//        Factory function to create rotational state propagator settings with generic stopping conditions.
//
//        Factory function to create rotational state propagator settings with generic stopping conditions.
//        It works by providing a key-value torque container, containing the list of torques acting on
//        each body. The map has as key a string denoting the name of the body on which a set of torques, provided
//        as value, act. This set is again a key-value container, with the key denoting the body exerting the torque
//        and the value (a pointer to) a torque settings object. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        torque_settings : SelectedTorqueMap
//            Set of torques acting on the bodies to propagate, provided as torque settings object.
//        bodies_to_integrate : list[str]
//            List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
//        initial_states : numpy.ndarray
//            Initial states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        propagator : RotationalPropagatorType, default=quaternions
//            Type of rotational propagator to be used (see `RotationalPropagatorType` enum).
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        RotationalStatePropagatorSettings
//            Rotational state propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "mass" && variant==0) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==1) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==2) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_settings : SelectedMassRateModelMap
//            Mass rates associated to each body, provided as a mass rate settings object.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==3) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==4) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
//        each body. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_settings : SelectedMassRateModelMap
//            Mass rates associated to each body, provided as a mass rate settings object.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "mass" && variant==0) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==1) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==2) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_settings : SelectedMassRateModelMap
//            Mass rates associated to each body, provided as a mass rate settings object.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==3) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==4) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
//        each body. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_settings : SelectedMassRateModelMap
//            Mass rates associated to each body, provided as a mass rate settings object.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "mass" && variant==0) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==1) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==2) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_settings : SelectedMassRateModelMap
//            Mass rates associated to each body, provided as a mass rate settings object.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==3) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==4) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
//        each body. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_settings : SelectedMassRateModelMap
//            Mass rates associated to each body, provided as a mass rate settings object.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "mass" && variant==0) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==1) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==2) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_settings : SelectedMassRateModelMap
//            Mass rates associated to each body, provided as a mass rate settings object.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==3) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==4) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
//        each body. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_settings : SelectedMassRateModelMap
//            Mass rates associated to each body, provided as a mass rate settings object.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "mass" && variant==0) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==1) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==2) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
//        each body. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_settings : SelectedMassRateModelMap
//            Mass rates associated to each body, provided as a mass rate settings object.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==3) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate models associated to
//        each body. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_models : dict[str, MassRateModel]
//            List of mass rates associated to each body, provided as mass rate models.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "mass" && variant==4) {
//            return R"(
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//
//        Factory function to create mass propagator settings with generic stopping conditions.
//        It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
//        each body. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        bodies_with_mass_to_propagate : list[str]
//            List of bodies whose mass should be numerically propagated.
//        mass_rate_settings : SelectedMassRateModelMap
//            Mass rates associated to each body, provided as a mass rate settings object.
//        initial_body_masses : numpy.ndarray
//            Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "multitype" && variant==0) {
//            return R"(
//
//        Factory function to create multitype propagator settings.
//
//        Factory function to create multitype propagator settings with generic stopping conditions.
//        It works by providing a list of SingleArcPropagatorSettings objects. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        propagator_settings_list : list[SingleArcPropagatorSettings]
//            List of SingleArcPropagatorSettings objects to use.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "multitype" && variant==1) {
//            return R"(
//
//        Factory function to create multitype propagator settings.
//
//        Factory function to create multitype propagator settings with generic stopping conditions.
//        It works by providing a list of SingleArcPropagatorSettings objects. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        propagator_settings_list : list[SingleArcPropagatorSettings]
//            List of SingleArcPropagatorSettings objects to use.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "multitype" && variant==0) {
//            return R"(
//
//        Factory function to create multitype propagator settings.
//
//        Factory function to create multitype propagator settings with generic stopping conditions.
//        It works by providing a list of SingleArcPropagatorSettings objects. In this function, the dependent variables to save are
//        provided as a unique DependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        propagator_settings_list : list[SingleArcPropagatorSettings]
//            List of SingleArcPropagatorSettings objects to use.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : DependentVariableSaveSettings, default=none
//            Cumulative dependent variable to be saved settings object (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//    } else if(name == "multitype" && variant==1) {
//            return R"(
//
//        Factory function to create multitype propagator settings.
//
//        Factory function to create multitype propagator settings with generic stopping conditions.
//        It works by providing a list of SingleArcPropagatorSettings objects. In this function, the dependent variables to save are provided
//        as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
//        through the termination settings object provided.
//
//
//        Parameters
//        ----------
//        propagator_settings_list : list[SingleArcPropagatorSettings]
//            List of SingleArcPropagatorSettings objects to use.
//        termination_settings : PropagationTerminationSettings
//            Generic termination settings object to check whether the propagation should be ended.
//        output_variables : list[SingleDependentVariableSaveSettings], default=[]
//            List of dependent variables to be saved (by default, no dependent variables are saved).
//        print_interval : float, default=TUDAT_NAN
//            Variable indicating how often (in seconds or in the unit of the independent variable) the current state and time are to be printed to the console (by default, they are never printed).
//
//        Returns
//        -------
//        MassPropagatorSettings
//            Mass propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "multi_arc" && variant==0) {
//            return R"(
//
//        Factory function to create multi-arc propagator settings.
//
//        Factory function to create multi-arc propagator settings. It works by providing separate settings for
//        each arc in a list.
//
//
//        Parameters
//        ----------
//        single_arc_settings : list[SingleArcPropagatorSettings]
//            List of SingleArcPropagatorSettings objects to use, one for each arc.
//        transfer_state_to_next_arc : bool, default=False
//            It denotes whether whether the initial state of arc N+1 is to be taken from arc N (for N>0).
//
//        Returns
//        -------
//        MultiArcPropagatorSettings
//            Multi-arc propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "hybrid_arc" && variant==0) {
//            return R"(
//
//        Factory function to create hybrid-arc propagator settings.
//
//        Factory function to create hybrid-arc propagator settings (i.e., a combination of single- and multi-arc dynamics).
//
//
//        Parameters
//        ----------
//        single_arc_settings : SingleArcPropagatorSettings
//            SingleArcPropagatorSettings object to use for the propagation.
//        multi_arc_settings : MultiArcPropagatorSettings
//            MultiArcPropagatorSettings object to use for the propagation.
//
//        Returns
//        -------
//        HybridArcPropagatorSettings
//            Hybrid-arc propagator settings object.
//
//    )";
//
//
//
//    } else if(name == "time_termination" && variant==0) {
//            return R"(
//
//        Factory function to create time termination settings for the propagation.
//
//        Factory function to create time termination settings for the propagation.
//        The propagation is stopped when the final time provided is reached.
//
//
//        Parameters
//        ----------
//        termination_time : float
//            Final time of the propagation.
//        terminate_exactly_on_final_condition : bool, default=False
//            Denotes whether the propagation is to terminate exactly on the final condition, or whether it is to terminate on the first step where it is violated.
//
//        Returns
//        -------
//        PropagationTimeTerminationSettings
//            Time termination settings object.
//
//    )";
//
//
//
//    } else if(name == "cpu_time_termination" && variant==0) {
//            return R"(
//
//        Factory function to create CPU time termination settings for the propagation.
//
//        Factory function to create CPU time termination settings for the propagation.
//        The propagation is stopped when the final CPU time provided is reached.
//
//
//        Parameters
//        ----------
//        cpu_termination_time : float
//            Maximum CPU time for the propagation.
//
//        Returns
//        -------
//        PropagationCPUTimeTerminationSettings
//            CPU time termination settings object.
//
//    )";
//
//
//
//    } else if(name == "dependent_variable_termination" && variant==0) {
//            return R"(
//
//        Factory function to create CPU time termination settings for the propagation.
//
//        Factory function to create CPU time termination settings for the propagation.
//        The propagation is stopped when the final CPU time provided is reached.
//
//
//        Parameters
//        ----------
//        dependent_variable_settings : SingleDependentVariableSaveSettings
//            Dependent variable object to be used as termination setting.
//        limit_value : float
//            Limit value of the dependent variable; if reached, the propagation is stopped.
//        use_as_lower_limit : bool, default=False
//            Denotes whether the limit value should be used as lower or upper limit.
//        terminate_exactly_on_final_condition : bool, default=False
//            Denotes whether the propagation is to terminate exactly on the final condition, or whether it is to terminate on the first step where it is violated.
//        termination_root_finder_settings : bool, default=None
//            Settings object to create root finder used to converge on exact final condition.
//
//        Returns
//        -------
//        PropagationDependentVariableTerminationSettings
//            Dependent variable termination settings object.
//
//    )";
//
//
//
//    } else if(name == "custom_termination" && variant==0) {
//            return R"(
//
//        Factory function to create custom termination settings for the propagation.
//
//        Factory function to create custom termination settings for the propagation.
//        The propagation is stopped when the condition provided is verified.
//
//
//        Parameters
//        ----------
//        custom_condition : Callable[[float], bool]
//            Function of time (independent variable) which is called during the propagation and returns a boolean value denoting whether the termination condition is verified.
//
//        Returns
//        -------
//        PropagationCustomTerminationSettings
//            Custom termination settings object.
//
//    )";
//
//
//
//    } else if(name == "hybrid_termination" && variant==0) {
//            return R"(
//
//        Factory function to create bybrid termination settings for the propagation.
//
//        Factory function to create hybrid termination settings for the propagation. This function can be used
//        to define that all conditions or a single condition of the conditions provided must be met to
//        stop the propagation.
//
//
//        Parameters
//        ----------
//        termination_settings : list[PropagationTerminationSettings]
//            List of single PropagationTerminationSettings objects to be checked during the propagation.
//        fulfill_single_condition : bool, default=False
//            Whether only a single condition of those provided must be met to stop the propagation (true) or all of them simultaneously (false).
//
//        Returns
//        -------
//        PropagationHybridTerminationSettings
//            Hybrid termination settings object.
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//
//
//
//
//namespace torque {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//    } else if(name == "AvailableTorque") {
//         return R"(
//
//        Enumeration of available torque types.
//
//        Enumeration of torque types supported by tudat.
//
//     )";
//
//
//
//    } else if(name == "TorqueSettings") {
//         return R"(
//
//        Functional base class to define settings for torques.
//
//     )";
//
//
//
//    } else if(name == "SphericalHarmonicTorqueSettings") {
//         return R"(
//
//        `TorqueSettings`-derived class to define settings for torques caused by spherical harmonic gravity.
//
//     )";
//
//
//
//
//    } else if(name == "aerodynamic" && variant==0) {
//            return R"(
//
//        Creates the settings for the aerodynamic torque.
//
//        Creates the settings for the aerodynamic torque exerted by a body with an atmosphere model and shape model on
//        another body. The body exerting the torque needs to have both an atmosphere model and a shape model defined.
//        Furthermore, the body undergoing the torque needs to have the aerodynamic coefficient interface and its moment
//        coefficients defined. In the case that the aerodynamic coefficients are defined as a function of the vehicle
//        orientation (e.g. angle of attack and sideslip angle), these angles can be manually or automatically defined.
//
//
//        Returns
//        -------
//        TorqueSettings
//            Torque settings object.
//
//    )";
//
//
//
//    } else if(name == "spherical_harmonic_gravitational" && variant==0) {
//            return R"(
//
//        Creates the settings for the spherical harmonic torque.
//
//        Torque exerted by a point mass on a body with an arbitrary degree/order spherical harmonics mass distribution.
//        The body exerting the torque only needs to have a gravitational model defined (point-mass or spherical harmonic),
//        while the body undergoing the torque needs to have a spherical harmonic gravity field defined.
//
//
//        Parameters
//        ----------
//        maximum_degree : int
//            Maximum degree of the spherical harmonic expansion.
//        maximum_order : int
//            Maximum order of the spherical harmonic expansion.
//
//        Returns
//        -------
//        TorqueSettings
//            Torque settings object.
//
//    )";
//
//
//
//    } else if(name == "second_degree_gravitational" && variant==0) {
//            return R"(
//
//        Creates the settings for the second-degree gravitational torque.
//
//        Torque exerted by a point mass on a body with a degree two spherical harmonics mass distribution.
//        A degree two spherical harmonics mass distribution can be represented by an inertia tensor; thus,
//        for this torque model, the body undergoing the torque needs to have an inertia tensor defined.
//        The body exerting the torque only needs to have a gravitational model defined (either point-mass or spherical
//        harmonics).
//
//
//        Returns
//        -------
//        TorqueSettings
//            Torque settings object.
//
//    )";
//
//
//
//    } else if(name == "custom" && variant==0) {
//            return R"(
//
//        Creates the settings for a custom torque.
//
//        Torque exerted by a point mass on a body with a degree two spherical harmonics mass distribution.
//        A degree two spherical harmonics mass distribution can be represented by an inertia tensor; thus,
//        for this torque model, the body undergoing the torque needs to have an inertia tensor defined.
//        The body exerting the torque only needs to have a gravitational model defined (either point-mass or spherical
//        harmonics).
//
//
//        Parameters
//        ----------
//        torque_function : Callable[[float], list]
//            Custom torque function with time as an independent variable.
//        scaling_function : Callable[[float], float], default=None
//            Scaling function with time as an independent variable to be multiplied by the custom torque function.
//
//        Returns
//        -------
//        TorqueSettings
//            Torque settings object.
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//
//
//
//}
//
//
//
//
//}
//
//
//
//
//
//namespace plotting {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//
//    } else if(name == "plot_blue_marble_ground_track" && variant==0) {
//            return R"(
//
//    )";
//
//
//
//    } else if(name == "plot_miller_ground_track" && variant==0) {
//            return R"(
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//
//
//
//
//namespace util {
//
//static inline std::string get_docstring(std::string name, int variant=0) {
//
//    if (name == "test") {
//        return "test";
//
//
//
//
//    } else if(name == "result2array" && variant==0) {
//            return R"(
//
//    )";
//
//
//
//    } else if(name == "compare_results" && variant==0) {
//            return R"(
//
//    )";
//
//
//
//    } else if(name == "redirect_std" && variant==0) {
//            return R"(
//
//    )";
//
//
//
//    } else {
//        return "No documentation found.";
//    }
//
//}
//
//
//}
//
//


}

