import { Simulation } from "./simulation.js";
import { TypeAnnee, c, k, h, G, AU, parsec, k_parsec, M_parsec, ly } from "./../../constants.js";
/**
 * @class Simulation_universe.
 * inheritance from Simulation class
 *
 * attributes :
 * @param temperature : current temperature of the universe.
 * @param hubble_cst : current value of the Hubble-LeMaître constant.
 * @param H0parsec : Hubble-Lemaître constant in international system units
 * @param matter_parameter : current value of the matter density parameter.
 * @param dark_energy : object containing current value of dark energy density parameter, value of w_0 and value of w_1.\
 * 	Note : When w_0 = -1 and w_0 = 0, the universe is equivalent to his counterpart with only a cosmologic constant.
 * @param constants : contains the value of the physics constants defined for the universe.
 * @param has_cmb : Has Cosmic Microwave Background (CMB).
 * @param has_neutrino : self explanatory.
 * @param is_flat : Forcing the curvature density parameter to 0.
 * @param is_single_matter : we use the single-fluid model with the matter dentity parameter equals to 1
 * @param is_single_cosmo : we use the single-fluid model with the dark energy density parameter equals to 1
 * @param is_single_radiation : we use the single-fluid model with the radiation dentity parameter equals to 1
 * @param is_single_radiation : we use the single-fluid model with the curvature dentity parameter equals to 1


 *
 * methods names :
 * @method modify_dark_energy
 * @method modify_constants
 * @method meter_to_light_year
 * @method meter_to_parsec
 * @method parsec_to_meters
 * @method seconds_to_years
 * @method runge_kutta_universe_1
 * @method runge_kutta_universe_2
 * @method calcul_rho_r
 * @method calcul_rho_lambda
 * @method calcul_rho_m
 * @method calcul_omega_r
 * @method calcul_omega_k
 * @method check_sum_omegas
 * @method Y
 * @method dY
 * @method F
 * @method function_E
 * @method T(z)
 * @method H(z)
 * @method omega_m_shift
 * @method omega_k_shift
 * @method omega_DE_shift
 * @method omega_r_shift
 * @method compute_scale_factor
 * @method compute_omegas
 * @method time
 * @method universe_age
 * @method duration
 * @method delta_dm
 * @method metric_distance
 * @method theta
 * @method theta_kpc
 * @method D
 * @method luminosity
 * @method luminosity_distance
 * @method light_distance
 * @method angular_diameter_distance
 * @method brightness
 * @method apparent_diameter
 * @method integral_duration_substituated
 * @method integral_distance
 * @method equa_diff_a
 * @method equa_diff_time
 */
export class Simulation_universe extends Simulation {
    //-------------------------constructor-----------------------
    constructor(id, temperature = 2.7255, hubble_cst = 67.74, matter_parameter = 3.089e-1, has_cmb = true, has_neutrino = true, is_flat = false) {
        super(id);
        this._dark_energy = {
            parameter_value: 6.911e-1,
            w_0: -1,
            w_1: 0,
        };
        this._constants = {
            c: c,
            k: k,
            h: h,
            G: G,
            AU: AU,
            parsec: parsec,
            k_parsec: k_parsec,
            M_parsec: M_parsec,
            ly: ly,
            nbrJours: TypeAnnee.Gregorienne,
        };
        this._temperature = temperature;
        this._hubble_cst = hubble_cst;
        this._H0parsec = (hubble_cst * 1e3) / (((AU * (180 * 3600)) / Math.PI) * 1e6);
        this._matter_parameter = matter_parameter;
        this._has_cmb = has_cmb;
        this._has_neutrino = has_neutrino;
        this._is_flat = is_flat;
        this.is_single_cosmo = false;
        this.is_single_curvature = false;
        this.is_single_matter = false;
        this.is_single_radiation = false;
    }
    //--------------------------Accessors------------------------
    // temperature
    get temperature() {
        return this._temperature;
    }
    set temperature(temperature) {
        this._temperature = temperature;
        this.check_sum_omegas();
    }
    // hubble_cst
    get hubble_cst() {
        return this._hubble_cst;
    }
    set hubble_cst(hubble_cst) {
        this._hubble_cst = hubble_cst;
        this._H0parsec=(hubble_cst * 1e3) / (((AU * (180 * 3600)) / Math.PI) * 1e6);
        this.check_sum_omegas();
    }
    //H0parsec
    get H0parsec() {
        return this._H0parsec;
    }
    set H0parsec(H0) {
        this._H0parsec = H0;
        this.check_sum_omegas();
    }
    // matter_parameter
    get matter_parameter() {
        return this._matter_parameter;
    }
    set matter_parameter(matter_parameter) {
        this._matter_parameter = matter_parameter;
        this.check_sum_omegas();
    }
    // dark_energy
    get dark_energy() {
        return this._dark_energy;
    }
    // constants
    get constants() {
        return this._constants;
    }
    // has_cmb
    get has_cmb() {
        return this._has_cmb;
    }
    set has_cmb(has_cmb) {
        this._has_cmb = has_cmb;
        this.check_sum_omegas();
    }
    // has_neutrino
    get has_neutrino() {
        return this._has_neutrino;
    }
    set has_neutrino(has_neutrino) {
        this._has_neutrino = has_neutrino;
        this.check_sum_omegas();
    }
    // is_flat
    get is_flat() {
        return this._is_flat;
    }
    set is_flat(is_flat) {
        this._is_flat = is_flat;
        this.check_sum_omegas();
    }
    get _is_single_matter() {
        return this._is_single_matter;
    }
    get _is_single_cosmo() {
        return this._is_single_cosmo;
    }
    get _is_single_curvature() {
        return this._is_single_curvature;
    }
    get _is_single_radiation() {
        return this._is_single_radiation;
    }
    //modifies the simulation if we use a single fluid model. Replaces setters for is_single_matter,
    //is_single_cosmo_cst, is_single_radiation and is_single_curvature
    single_fluid(model) {
        let w0 = this.dark_energy.w_0;
        let w1 = this.dark_energy.w_1;
        let T = 0;
        let omega_m = 0;
        let omega_lambda = 0;
        this.is_single_cosmo=false;
        this.is_single_curvature=false;
        this.is_single_matter=false;
        this.is_single_radiation=false;
        if (model == "matter") {
            omega_m = 1;
            this.is_single_matter = true;
        }
        if (model == "cosmo_cst") {
            omega_lambda = 1;
            this.is_single_cosmo = true;
        }
        if (model == "radiation") {
            this.is_single_radiation = true;
            omega_r = 1;
            let c = this.constants.c;
            let h = this.constants.h;
            let pi = Math.PI;
            let G = this.constants.G;
            let k = this.constants.k;
            let H = this._H0parsec;
            T = Math.pow(45 * Math.pow(c, 5) * Math.pow(h, 3) / (64 * Math.pow(pi, 6) * G * Math.pow(k, 4)), 1 / 4) * Math.pow(H, 1 / 2);
            let A = 45 * Math.pow(c, 2) * Math.pow(h, 3);
            let B = 64 * Math.pow(pi, 6) * G * Math.pow(k, 4);
            T = Math.pow(A / B, 1 / 4) * Math.pow(H, 1 / 2);
        }
        if (model == "curvature") {
            this.is_single_curvature = true;
        }
        this.matter_parameter = omega_m;
        this.modify_dark_energy(omega_lambda, w0, w1);
        this.calcul_omega_r();
        this.temperature = T;
        this.calcul_omega_k();
        this.check_sum_omegas();
    }
    //---------------------------methods-------------------------
    //                      redefined methods
    //                         new methods
    /**
     * replace the setter for the dark_energy attribute
     *
     * @param DE_parameter_value value of the dark energy density parameter
     * @param DE_w_0 value of w_0
     * @param DE_w_1 value of w_1
     *
     * Note : w_0, w_1 are parameters that describe the nature of the dark energy.
     */
    modify_dark_energy(DE_parameter_value, DE_w_0, DE_w_1) {
        if (DE_parameter_value !== undefined) {
            this._dark_energy.parameter_value = DE_parameter_value;
            this.check_sum_omegas(true);
        }
        if (DE_w_0 !== undefined) {
            this._dark_energy.w_0 = DE_w_0;
        }
        if (DE_w_1 !== undefined) {
            this._dark_energy.w_1 = DE_w_1;
        }
    }
    /**
     * replace the setter for the constants attribute
     *
     * @param c light speed constant
     * @param k Boltzmann constant
     * @param h Planck constant
     * @param G Newton constant
     * @param TypeAnnee Number of days in chosen Type of Year
     */
    modify_constants(c, k, h, G, typeAnnee) {
        if (c !== undefined) {
            this._constants.c = c;
        }
        if (k !== undefined) {
            this._constants.k = k;
        }
        if (h !== undefined) {
            this._constants.h = h;
        }
        if (G !== undefined) {
            this._constants.G = G;
        }
        if (typeAnnee = "Sidérale") {
            var nbrjours = TypeAnnee.Siderale;
        }
        else if (typeAnnee = "Julienne") {
            var nbrjours = TypeAnnee.Julienne;
        }
        else if (typeAnnee = "Tropique (2000)") {
            var nbrjours = TypeAnnee.Tropique2000;
        }
        else {
            var nbrjours = TypeAnnee.Gregorienne;
        }
    }
    /**
     * Converts a length in meters to a lenght in light years
     * @param l : length in meters
     */
    meter_to_light_year(l) {
        return Number(l) / this.constants.ly;
    }
    /**
 * Converts a length in meters to a lenght in parsec
 * @param l : length in meters
 */
    meter_to_parsec(l) {
        return Number(l) / this.constants.parsec;
    }
    /**
     * Converts a length in meters to a lenght in kiloparsec
     * @param l : length in meters
     */
    meter_to_kiloparsec(l) {
        return Number(l) / this.constants.k_parsec;
    }
    /**
*Converts a length in parsec to a lenght in meters
@param l : length in parsec
    */
    parsec_to_meters(l) {
        return Number(l) * this.constants.parsec;
    }
    /**
*Converts a length in kiloparsec to a lenght in meters
@param l : length in parsec
    */
    kiloparsec_to_meters(l) {
        return Number(l) * this.constants.k_parsec;
    }
    /**
     * Convert a duration in seconds to a duration in years
     */
    seconds_to_years(t) {
        return Number(t) / (Number(this.constants.nbrJours) * 24 * 60 * 60);
    }
    /**
     * Fourth order Runge-Kutta method for second order derivatives for universe computation.
     *
     * @param step The step of computation
     * @param x_0 x_point where the calcul start
     * @param y_0 initial value of y at x_0
     * @param interval Array containing [xmin, xmax]
     * @param funct function or method that define the equation to resolve, your function has to accept 2 numbers and return a number
     *
     * @returns [step: number, x: number[], y:number[]].
     */
    runge_kutta_universe_1(step, x_0 = 0, y_0 = 1, funct, interval = [0, 5]) {
        // Init parameter
        let x = [x_0];
        let y = [y_0];
        // Computation loops
        // Computing with a positive step, i increments the array
        let i = 0;
        let result_runge_kutta;
        while (interval[0] <= x[i] && x[i] < interval[1]) {
            result_runge_kutta = this.runge_kutta_equation_order1(this, step, x[i], y[i], funct);
            x.push(result_runge_kutta[0]);
            y.push(result_runge_kutta[1]);
            i++;
        }
        /*
            Computing with a negative step,
            since we decrease the value of x we add the elements at the beginning of the arrays,
            so for each step we take the first element of the array to compute the next one.
        */
        while (interval[0] <= x[0] && x[0] < interval[1]) {
            result_runge_kutta = this.runge_kutta_equation_order1(this, -step, x[0], y[0], funct);
            x.unshift(result_runge_kutta[0]);
            y.unshift(result_runge_kutta[1]);
            i++;
        }
        return {
            x: x,
            y: y
        };
    }
    /**
     * Fourth order Runge-Kutta method for second order derivatives for universe computation.
     *
     * @param step The step of computation
     * @param x_0 x_point where the calcul start
     * @param y_0 initial value of y at x_0
     * @param yp_0 initial value of the derivative of y at x_0
     * @param interval Array containing [ymin, ymax]
     * @param funct function or method that define the equation to resolve, your function has to accept 3 numbers and return a number
     *
     * @returns [step: number, x: number[], y:number[], yp: number[]].
     */
    runge_kutta_universe_2(step, x_0 = 0, y_0 = 1, dy_0 = 1, funct, interval) {
        // Init parameter
        let x = [x_0];
        let y = [y_0];
        let dy = [dy_0];
        // Computation loops
        // Computing with a positive step, i increments the array
        let i = 0;
        let result_runge_kutta;
        while (interval[0] <= y[i] && y[i] < interval[1]) {
            result_runge_kutta = this.runge_kutta_equation_order2(this, step, x[i], y[i], dy[i], funct);
            x.push(result_runge_kutta[0]);
            y.push(result_runge_kutta[1]);
            dy.push(result_runge_kutta[2]);
            i++;
        }
        /*
            Computing with a negative step,
            since we decrease the value of x we add the elements at the beginning of the arrays,
            so for each step we take the first element of the array to compute the next one.
        */
        while (interval[0] <= y[0] && y[0] < interval[1]) {
            result_runge_kutta = this.runge_kutta_equation_order2(this, -step, x[0], y[0], dy[0], funct);
            x.unshift(result_runge_kutta[0]);
            y.unshift(result_runge_kutta[1]);
            dy.unshift(result_runge_kutta[2]);
            i++;
        }
        return {
            x: x,
            y: y,
            dy: dy
        };
    }
    calcul_rho_r() {
        let sigma = (2 * Math.pow(Math.PI, 5) * Math.pow(this.constants.k, 4)) /
            (15 * Math.pow(this.constants.h, 3) * Math.pow(this.constants.c, 2));
        return (4 * sigma * Math.pow(this.temperature, 4)) / Math.pow(this.constants.c, 3);
    }
    calcul_rho_lambda() {
        let omega = this.dark_energy.parameter_value;
        let const_cosmo = 3 * Math.pow(this._H0parsec, 2) * omega / Math.pow(this.constants.c, 2);
        return const_cosmo * Math.pow(c, 2) / (8 * Math.PI * G);
    }
    calcul_rho_m() {
        return 3 * Math.pow(this.hubble_cst * 3.086 * Math.pow(10, -16), 2) / (8 * Math.PI * G);
    }
    /**
     * compute radiation density parameter at current time
     * @returns the radiation density parameter
     *
     * 		if (universe_age === undefined) {
            age = this.universe_age();
     */
    calcul_omega_r() {
        // Hubble-Lemaître constant in international system units (Système International)
        let omega_r = (8 * Math.PI * this.constants.G * this.calcul_rho_r()) / (3 * Math.pow(this._H0parsec, 2));
        if (this.has_neutrino) {
            omega_r *= 1.68;
        }
        if (this.is_single_radiation) {
            omega_r = 1;
        }
        if (!this.has_cmb || this.is_single_cosmo || this.is_single_curvature || this.is_single_matter) {
            omega_r = 0;
        }
        return omega_r;
    }
    /**
     * Compute curvature density parameter at current time
     * @returns the curvature density parameter
     */
    calcul_omega_k() {
        if (this.is_flat) {
            return 0;
        }
        else {
            return (1 -
                this.calcul_omega_r() -
                this.matter_parameter -
                this.dark_energy.parameter_value);
        }
    }
    /**
     * Check if the sum of the density parameters is equal to 1. Otherwise modify one parameter to correct the sum.
     * @param modify_matter true : modify the matter parameter, false : dark energy parameter instead
     * @returns false if one parm has been modified, true otherwise
     */
    check_sum_omegas(modify_matter = true) {
        let is_param_modified = false;
        let omega_r = this.calcul_omega_r();
        let sum = this.matter_parameter + omega_r + this.dark_energy.parameter_value + this.calcul_omega_k();
        if (this.is_flat && sum !== 1) {
            is_param_modified = true;
            if (!modify_matter) {
                this.matter_parameter = 1 - this.dark_energy.parameter_value - omega_r;
            }
            else {
                this.modify_dark_energy(1 - this.matter_parameter - omega_r);
            }
        }
        return is_param_modified;
    }
    /**
     * Y function \
     * see Theory about cosmology and dark_energy
     * @param x variable
     * @returns value of Y at position x
     */
    Y(x) {
        return Math.exp(-3 *
            (this.dark_energy.w_0 + this.dark_energy.w_1 + 1) *
            Math.log(Number(x)) -
            3 * this.dark_energy.w_1 * (1 - Number(x)));
    }
    /**
     * Y' function \
     * see Theory about cosmology and dark_energy
     * @param x variable
     * @returns value of the derivative of Y at position x
     */
    dY(x) {
        return (this.Y(x) *
            (3 * this.dark_energy.w_1 -
                3 * (1 + this.dark_energy.w_0 + this.dark_energy.w_1) / x));
    }
    /**
     * F function \
     * see Theory about cosmology and dark_energy
     * @param x variable
     * @returns value of F(x)
     */
    F(x) {
        return (Math.pow((1 + x), 2) * this.calcul_omega_k() +
            Math.pow((1 + x), 3) * this.matter_parameter +
            Math.pow((1 + x), 4) * this.calcul_omega_r() +
            this.Y(1 / (1 + x)) * this.dark_energy.parameter_value);
    }
    /**
     * E function \
     * will be used tu calculate the Omegas in function of the shift z
     * @param x
     * @param omegam0 matter density parameter
     * @param omegalambda0 dark energy parameter
     * @param Or radiation density parameter
     * @returns
     */
    function_E(x, omegam0, omegalambda0, Or) {
        let w0 = this.dark_energy.w_0;
        let w1 = this.dark_energy.w_1;
        let Yde = Math.exp(-3 * (w1 + w0 + 1) * Math.log(1 / (1 + x)) - (3 * w1 * (1 - (1 / (1 + x)))));
        return (Number(Or) * Math.pow((1 + x), 4) + Number(omegam0) * Math.pow((1 + x), 3)
            + (1 - Number(omegam0) - Number(Or) - Number(omegalambda0)) * Math.pow((1 + x), 2)
            + Number(omegalambda0) * Yde);
    }
    /**
     * will be used for reverse calculations
     */
    derivative_function_E(x) {
        let w0 = this.dark_energy.w_0;
        let w1 = this.dark_energy.w_1;
        let omega_m = this.matter_parameter;
        let omega_r = this.calcul_omega_r();
        let omega_lambda = this.dark_energy.parameter_value;
        let U = 3 * (w1 / Math.pow(1 + Number(x), 2) + (w0 + w1 + 1) / (1 + Number(x)));
        return U * this.function_E(Number(x), omega_m, omega_lambda, omega_r);
    }
    /**
     * calculates the temperature as a function of the shift z
     * @param z shift
     *
     */
    T(z) {
        return this.temperature * (1 + Number(z));
    }
    /**
     * calculates the Hubble constant as a function of the shift z
     * @param z
     */
    H(z) {
        let omega_m0 = this.matter_parameter;
        let omega_DE0 = this.dark_energy.parameter_value;
        let omega_r0 = this.calcul_omega_r();
        return this.hubble_cst * Math.pow(this.function_E(Number(z), omega_m0, Number(omega_DE0), omega_r0), 0.5);
    }
    /**
     * calculates the matter density parameter as a function of the shif z
     * @param z Redshift
     * */
    omega_m_shift(z) {
        //		return this.matter_parameter * (1 + z)**3 / this.F(z);
        let omega_m0 = this.matter_parameter;
        let omega_DE0 = this.dark_energy.parameter_value;
        let omega_r0 = this.calcul_omega_r();
        return omega_m0 * Math.pow(1 + Number(z), 3) / this.function_E(Number(z), omega_m0, Number(omega_DE0), omega_r0);
    }
    /**
 * calculates the curvature density parameter as a function of the shif z
 * @param z Redshift
 * */
    omega_k_shift(z) {
        let omega_k0 = this.calcul_omega_k();
        let omega_m0 = this.matter_parameter;
        let omega_DE0 = this.dark_energy.parameter_value;
        let omega_r0 = this.calcul_omega_r();
        return omega_k0 * Math.pow(1 + Number(z), 2) /
            this.function_E(Number(z), omega_m0, Number(omega_DE0), omega_r0);
    }
    /**
    * calculates the dark energy parameter as a function of the shif z
    * @param z Redshift
    * */
    omega_DE_shift(z) {
        let omegaDE0 = this.dark_energy.parameter_value;
        let omega_m0 = this.matter_parameter;
        let omega_r0 = this.calcul_omega_r();
        return omegaDE0 / this.function_E(Number(z), omega_m0, omegaDE0, omega_r0);
    }
    /**
    * calculates the radiation density parameter as a function of the shif z
    * @param z Redshift
    * */
    omega_r_shift(z) {
        let omega_r0 = this.calcul_omega_r();
        let omega_m0 = this.matter_parameter;
        let omega_DE0 = this.dark_energy.parameter_value;
        return omega_r0 * Math.pow(1 + Number(z), 4) / this.function_E(Number(z), omega_m0, omega_DE0, omega_r0);
        //return this.calcul_omega_r() * (1 + z)**3 / this.F(z);
    }
    /**
     * compute the scale factor of the universe as function of time
     * @param step Computation step
     * @param interval_a Array containing a_min et a_max value
     * @param universe_age Permit to pass an already computed value for the universe age. If not given, the method recompute the value.
     * @returns t value, a value, derivative of a
     */
    compute_scale_factor(step, interval_a = [0, 5], universe_age) {
        let age;
        if (universe_age === undefined) {
            age = this.universe_age();
        }
        else {
            age = universe_age;
        }
        if (isNaN(age)) {
            age = 0;
        }
        let result = this.runge_kutta_universe_2(step, 0, 1, 1, this.equa_diff_a, interval_a);
        for (let index = 0; index < result.x.length; index++) {
            result.x[index] = (result.x[index] / this._H0parsec + age) / (3600 * 24 * 365.2425);
        }
        return result;
    }
    /**
     * Computing the 4 density parameters given an array of cosmologic shift value
     * @param z_array array containing z points where to compute the omegas
     */
    compute_omegas(z_array) {
        let omega_matter = [];
        let omega_rad = [];
        let omega_de = [];
        let omega_courbure = [];
        let radiation = this.calcul_omega_r();
        let curvature = this.calcul_omega_k();
        z_array.forEach(z => {
            omega_matter.push(this.matter_parameter * Math.pow((1 + z), 3) / this.F(z));
            omega_rad.push(radiation * Math.pow((1 + z), 3) / this.F(z));
            omega_de.push(this.dark_energy.parameter_value * Math.pow((1 + z), 3) / this.F(z));
            omega_courbure.push(curvature * Math.pow((1 + z), 3) / this.F(z));
        });
        return {
            omega_matter: omega_matter,
            omega_rad: omega_rad,
            omega_de: omega_de,
            omega_courbure: omega_courbure
        };
    }
    /**
     * Compute the time as a function of the cosmologic shift
     * @param n number of computation points
     * @param zmin
     * @param zmax
     * @returns
     */
    time(n, zmin, zmax) {
        let step = (zmax - zmin) / n;
        let time_zmin;
        try {
            time_zmin = this.duration(0, zmin);
        }
        catch (e) {
            return e;
        }
        return this.runge_kutta_universe_1(step, zmin, time_zmin, this.equa_diff_time, [zmin, zmax]);
    }
    universe_age() {
        /*
        To compute the age of the universe we need to integrate from x = 0 to x -> infinity. To resolve this problem we do a substitution with
        x = y / (1 - y) which implies dx = dy / (1 - y)². This result with an integral from y = 0 to y = 1 that can be digitally resolved.
        */
        let age = 1 / this._H0parsec;
        if (this.is_single_matter) {
            age *= 2 / 3;
        }
        else if (this.is_single_radiation) {
            age *= 1 / 2;
        }
        else if (this.is_single_cosmo) {
            return NaN;
        }
        else if (!this.is_single_radiation) {
            age *=
                this.simpson(this, this.integral_duration_substituated, 0, 1, 10000);
        }
        return age;
    }
    /**
     * name
     */
    emission_age(z) {
        let infimum = z / (1 + z);
        let age;
        let H = this._H0parsec;
        if (this.is_single_cosmo) {
            age = NaN;
        }
        else if (this.is_single_curvature) {
            age = 1 / (H * (1 + z));
        }
        else if (this.is_single_matter) {
            age = (2 / 3) * Math.pow(1 + z, -3 / 2) / H;
        }
        else if (this.is_single_radiation) {
            age = (1 / 2) * (1 / H) * Math.pow(1 + z, -2);
        }
        else {
            age =
                this.simpson(this, this.integral_duration_substituated, infimum, 1, 1000) / H;
        }
        return age;
    }
    /**
     * Compute the cosmologic duration between two cosmologics shift z
     * @param z_1 the closest cosmologic shift from ours (z = 0)
     * @param z_2 the farest cosmologic shift from ours (z = 0)
     * @returns error if z_1 or z_2 < -1, duration if both value are accepted.
     */
    duration(z_1, z_2) {
        if (z_1 <= -1 || z_2 <= -1) {
            throw new Error("Cosmologic shift z cannot be equal or lower than -1 included");
        }
        let duration;
        if (this.is_single_curvature || this.is_single_matter || this.is_single_radiation) {
            duration = this.emission_age(z_2) - this.emission_age(z_1);
        }
        else if (this.is_single_cosmo) {
            duration = (1 / this.H0parsec) * Math.log(1 + z_1 / (1 + z_2));
        }
        let infimum = z_1 / (1 + z_1);
        let supremum = z_2 / (1 + z_2);
        duration = this.simpson(this, this.integral_duration_substituated, infimum, supremum, 1000) / this._H0parsec;
        return duration;
    }
    /**
     * Compute the distance between two cosmologics shift z
     * @param z_1 the closest cosmologic shift from ours (z = 0)
     * @param z_2 the farest cosmologic shift from ours (z = 0)
     * @returns error if z_1 or z_2 < -1, duration if both value are accepted.
     */
    delta_dm(z1, z2) {
        let distance;
        let curvature = this.calcul_omega_k();
        distance = this.simpson(this, this.integral_distance, Number(z1), Number(z2), 100000);
        if (curvature < 0) {
            distance =
                Math.sinh(Math.sqrt(Math.abs(curvature)) * distance) /
                    Math.sqrt(Math.abs(curvature));
        }
        else if (curvature > 0) {
            distance =
                Math.sin(Math.sqrt(Math.abs(curvature)) * distance) /
                    Math.sqrt(Math.abs(curvature));
        }
        distance *= this.constants.c / this._H0parsec;
        return distance;
    }
    /**
     * Compute the distance between us and an object at a cosmologic redshit z
     * @param z cosmologic shift
     * @returns the distance
     */
    metric_distance(z) {
        let distance = this.delta_dm(0, z);
        return distance;
    }
    /**
     * Computes the value of theta in second of arc
     * @param  D
     * @param  z Cosmologic shift
     * @param dm Metric distance
     * @returns theta
     */
    theta(D, z, dm) {
        let distance;
        if (dm === undefined) {
            distance = Number(this.metric_distance(z));
        }
        else {
            distance = Number(dm);
        }
        let add = this.angular_diameter_distance(z, distance);
        return 206265 * Number(D) / add;
    }
    /**
 * Computes the value of theta in second of arc if D is in kiloparsec
 * @param  D
 * @param  z Cosmologic shift
 * @param dm Metric distance
 * @returns theta
 */
    theta_kpc(D, z, dm) {
        let D_m = this.kiloparsec_to_meters((D));
        return this.theta(D_m, z, dm);
    }
    /**
    * Computes the value of D in meters and in parsec
    * @param  theta (in seconds of arcs)
    * @param  z Cosmologic shift
    * @param dm Metric distance
    * @returns D
    */
    D(theta, z, dm) {
        let distance;
        if (dm === undefined) {
            distance = Number(this.metric_distance(z));
        }
        else {
            distance = Number(dm);
        }
        let add = this.angular_diameter_distance(z, distance);
        let D_meters = add * Number(theta) / 206265;
        let D_kpc = this.meter_to_kiloparsec(D_meters);
        return { "m": D_meters, "kpc": D_kpc };
    }
    /**
     * Compute the luminosity distance
     * @param z Cosmologic shift
     * @param distance_metric optionnal parameters for optimisation (permit you to pass an already calculated distances for optimisation)
     * @returns luminosity distance
     */
    luminosity_distance(z, distance_metric) {
        let distance;
        if (distance_metric === undefined) {
            distance = Number(this.metric_distance(z));
        }
        else {
            distance = distance_metric;
        }
        let lum_dis = distance * (1 + z);
        // dictionnaire test 
        let values_unit = { pc: this.meter_to_parsec(lum_dis), meter: lum_dis, ly: this.meter_to_light_year(lum_dis) };
        return values_unit;
    }
    /**
     * @param z Cosmologic shift
     * @param distance_metric
     */
    angular_diameter_distance(z, distance_metric) {
        let distance;
        if (distance_metric === undefined) {
            distance = Number(this.metric_distance(z));
        }
        else {
            distance = distance_metric;
        }
        return distance / (1 + z);
    }
    /**
     * This fonction compute the distance with the simple formula c*t.
     * @param z Cosmologic shift
     * @returns c*t
     */
    light_distance(z) {
        let duration = this.duration(0, z);
        let c = this.constants.c;
        return duration * c;
    }
    /**
     * Compute the luminosity of an astronomical object of an unifrom intensity I
     * @param I intensity
     * @returns luminosity
     */
    luminosity(I) {
        return 4 * Math.PI * I;
    }
    /**
     * Compute the brightness of an object situated at a cosmologic redshit z
     * @param z Cosmologic shift
     * @param luminosity self explanatory
     * @param distance_metric optionnal parameters for optimisation (permit you to pass an already calculated distances for optimisation)
     */
    brightness(z, luminosity, distance_metric) {
        let distance;
        if (distance_metric === undefined) {
            distance = Number(this.metric_distance(z));
        }
        else {
            distance = distance_metric;
        }
        return luminosity / (4 * Math.PI * Math.pow(distance * (1 + z), 2));
    }
    /**
     * Compute the apparent diameter (Or the angle between 2 object of same shift z)
     * @param D_e Euclydien linear diameter
     * @param z Cosmologic shift
     * @param distance_metric optionnal parameters for optimisation (permit you to pass an already calculated distances for optimisation)
     * @returns The apparent diameter
     */
    apparent_diameter(D_e, z, distance_metric) {
        let distance;
        if (distance_metric === undefined) {
            distance = Number(this.metric_distance(z));
        }
        else {
            distance = distance_metric;
        }
        let values_unit2 = { meter: (D_e * (1 + z)) / distance, pc: this.meter_to_parsec((D_e * (1 + z)) / distance), ly: this.meter_to_light_year((D_e * (1 + z)) / distance) };
        return values_unit2;
    }
    /**
     * formula 1/(1 + x) * sqrt(1 / F)
     * @param Simu object in witch method is applied, permit to use this function with simulation method
     * @param x variable
     * @returns 1/(1 + x) * sqrt(1 / F)
     */
    integral_duration(Simu, x) {
        return ((1 / (1 + x)) * Math.sqrt(1 / Simu.F(x)));
    }
    /**
     * formula 1/(1 + x) * 1/sqrt(F) with the substitution x = y/(1 - y), to be used with simpson method to compute duration.
     * @param Simu object in witch method is applied, permit to use this function with simulation method
     * @param y variable
     * @returns (1 - y) * 1/sqrt(F(x)) * 1/(1 - y)²\
     *
     * Note : 1/(1 - y)² is the term come from dx = dy/(1 - y)²
     */
    integral_duration_substituated(Simu, y) {
        return ((((1 - y) / Math.pow(1 - y, 2))) /
            Math.sqrt(Simu.F(y / (1 - y))));
    }
    /**
     * Integral used to compute the distances
     * @param Simu object in witch method is applied, permit to use this function with simulation method
     * @param x variable
     * @returns 1/F²(x)
     */
    integral_distance(Simu, x) {
        return 1 / Math.sqrt(Simu.F(x));
    }
    /**
     * Right part of the differential equation of a(tau) designed to be used in runge_kutta_universe_2 method
     * @param Simu object in witch method is applied, permit to use this function with simulation method
     * @param tau time
     * @param a function a(t)
     * @param da derivative of a(t)
     * @returns result of the right part\
     * Note: tau and da are not used but have to be defined for this method to be accepted in the runge_kutta_equation_order2 method of simulation class
     */
    equa_diff_a(Simu, tau, a, da = 0) {
        let omega_r = Simu.calcul_omega_r();
        let omega_m = Simu.matter_parameter;
        let omega_de = Simu.dark_energy.parameter_value;
        return (-(omega_r / Math.pow(a, 3)) -
            (0.5 * omega_m) / Math.pow(a, 2) +
            omega_de *
                (a * Simu.Y(a) + (Math.pow(a, 2) * Simu.dY(a)) / 2));
    }
    /**
     * Right part of the differential equation of t(z) designed to be used in runge_kutta_universe_1 method
     * @param Simu object in witch method is applied, permit to use this function with simulation method
     * @param z Cosmologic shift
     * @param t function time t(z)
     * @returns result of the right part\
     * Note: t is not used but has to be defined for this method to be accepted in the runge_kutta_equation_order1 method of simulation class
     */
    equa_diff_time(Simu, z, t = 0) {
        return 1 / (this.H0parsec * (1 + z) * Math.sqrt(Simu.F(z)));
    }
    //CALCULS INVERSES
    //REVERSE CALCULATIONS 
    reverse_shift_dm(dm_param) {
        let dm = Number(dm_param);
        let omega_k = this.calcul_omega_k();
        let za = 0;
        let zb = 5e10;
        let ex = 0.00001; //indicates the permissible uncertainty for the dichotomie method
        //"reconditionner" is a list of two values : [new_zb,constraint]. constraint=0 or 1.
        // 0 means there is no constraint.
        let reconditionner = this.reconditionner(za, zb);
        zb = Number(reconditionner[0]) - 0.0001;
        let constraint = Number(reconditionner[1]);
        let dm_za = this.metric_distance(za);
        let dm_zb = this.metric_distance(zb);
        if (Number(dm) === 0) {
            return 0;
        }
        if (omega_k <= 0) {
            let limit = 0;
            if (constraint == 0) {
                while (dm > dm_zb && limit < 100) {
                    zb = zb * 10;
                    dm_zb = this.metric_distance(zb);
                    limit += 1;
                }
                if (limit >= 100) {
                    return NaN;
                }
            }
            /* the method "metric_distance" is monotonous if omega_k<=0. If dm > dm_zb and that
            zb is at a maximum, no solution can be found*/
            if (dm > dm_zb) {
                return NaN;
            }
            let Z = this.dichotomie(za, zb, this.metric_distance_simpson, dm, ex);
            return Z;
        }
        else {
            //Amplitude A of the method "delta_dm"
            let A = (this.constants.c / (this._H0parsec * Math.sqrt(Math.abs(omega_k))));
            if (dm > A) {
                return NaN;
            }
            let integB = Math.sqrt(Math.abs(omega_k)) * this.simpson(this, this.metric_distance_simpson, 0, zb, 10000);
            //this.simpson(0, zb, fonction_dm, omegam0, Number(omegalambda0), Number(Or), eps);
            //In this case, sin(integrale) is monotonous on the interval [za;zb]
            if (integB < Math.PI / 2) {
                //We check that dm<dm_zb because dm_zb<A.
                if (dm > dm_zb) {
                    return NaN;
                }
                return this.dichotomie(0, zb, this.metric_distance_simpson, dm, ex);
            }
            else if ((integB > Math.PI / 2) && (integB < Math.PI)) {
                let z_Pi_div_2 = this.dichotomie(0, zb, this.Integral_dm, Math.PI / 2, ex);
                let z_sol_1 = this.dichotomie(0, z_Pi_div_2, this.metric_distance_simpson, dm, ex);
                if (dm > dm_zb) {
                    this.dichotomie(z_Pi_div_2, zb, this.metric_distance_simpson, dm, ex);
                    return (z_sol_1 + ", " + z_sol_2);
                }
                return z_sol_1;
            }
            else {
                let z_Pi_div_2 = Number(dichotomie(0, zb, this.Integral_dm, Math.PI / 2, ex));
                let z_Pi = Number(this.dichotomie(0, 5e10, this.Integral_dm, Math.PI, ex));
                let z_sol_1 = this.dichotomie(0, z_Pi_div_2, this.metric_distance_simpson, dm, ex);
                let z_sol_2 = this.dichotomie(z_Pi_div_2, z_Pi, this.metric_distance_simpson, dm, ex);
                let z_f2 = z_sol_1 + ", " + z_sol_2;
                return z_f2;
            }
        }
    }
    reconditionner(za, zb) {
        let ex = 0.000001;
        let omega_r = this.calcul_omega_r();
        let omega_m = this.matter_parameter;
        let omega_lambda = this.dark_energy.parameter_value;
        let a = omega_r;
        let b = 4 * omega_r + omega_m;
        //We search the solutions of a polynomial of degree four
        if (a != 0 && (Math.abs(a / b) > 1e-3 || Math.abs(a / c) > 1e-3)) {
            let roots = this.fourth_order_solver();
            if (roots.length >= 2) {
                //We have at least two positive solutions. We select the smaller one.
                if ((roots[0] > 0) && (roots[1] > 0)) {
                    return [roots[1], 1];
                }
                //We hace at least two solution. One is positive and one is negative. We select the positive one.
                else if ((roots[0] > 0) && (roots[1] < 0)) {
                    return [roots[0], 1];
                }
            }
        }
        //We search the solutions of a polynomial of degree three
        else {
            let D_prime_za = this.derivative_function_E(za);
            let D_prime_zb = this.derivative_function_E(zb);
            if (D_prime_za * D_prime_zb < 0) {
                let z_prime = Number(this.dichotomie(0, zb, this.derivative_function_E, 0, ex));
                if (this.function_E(z_prime, omega_m, omega_lambda, omega_r) == 0) {
                    return [z_prime, 1];
                }
                else if (this.function_E(z_prime, omega_m, omega_lambda, omega_r) < 0) {
                    let z_pr = Number(this.dichotomie(0, z_prime, this.function_E, 0, ex));
                    return [z_pr, 1];
                }
            }
            else {
                if (this.function_E(za, omega_m, omega_lambda, omega_r) * this.function_E(zb, omega_m, omega_lambda, omega_r) < 0) {
                    let z_pr = Number(this.dichotomie(za, zb, this.function_E, 0, ex));
                    return [z_pr, 1];
                }
            }
        }
        return [zb, 0];
    }
    fourth_order_solver() {
        let omega_m = this._matter_parameter;
        let omega_lambda = this.dark_energy.parameter_value;
        let omega_r = this.calcul_omega_r();
        let a = omega_r;
        let b = 4 * omega_r + omega_m;
        let C = 5 * omega_r + 2 * omega_m - omega_lambda + 1;
        let d = 2 * omega_r + omega_m - 2 * omega_lambda + 2;
        let p = -3 * Math.pow(b / a, 2) / 8 + C / a;
        let q = Math.pow(b / a, 3) / 8 - (C / a) * (b / a) / 2 + d / a;
        let U1 = this.cubic_root_searcher();
        let z1 = (-Math.pow(U1 - p, 0.5) + Math.pow((U1 - p) - 4 * (0.5 * U1 - q / (2 * Math.pow(U1 - p, 0.5))), 0.5)) / 2;
        let z2 = (-Math.pow(U1 - p, 0.5) - Math.pow((U1 - p) - 4 * (0.5 * U1 - q / (2 * Math.pow(U1 - p, 0.5))), 0.5)) / 2;
        let z3 = (Math.pow(U1 - p, 0.5) + Math.pow((U1 - p) - 4 * (0.5 * U1 + q / (2 * Math.pow(U1 - p, 0.5))), 0.5)) / 2;
        let z4 = (Math.pow(U1 - p, 0.5) - Math.pow((U1 - p) - 4 * (0.5 * U1 + q / (2 * Math.pow(U1 - p, 0.5))), 0.5)) / 2;
        let x1 = z1 - 0.25 * (b / a);
        let x2 = z2 - 0.25 * (b / a);
        let x3 = z3 - 0.25 * (b / a);
        let x4 = z4 - 0.25 * (b / a);
        let array_roots = [x1, x2, x3, x4];
        for (var i = array_roots.length - 1; i >= 0; i--) {
            if (isNaN(array_roots[i])) {
                array_roots.splice(i, 1);
            }
        }
        array_roots.sort(function (a, b) { return b - a; });
        return array_roots;
    }
    cubic_root_searcher() {
        let omega_m = this._matter_parameter;
        let omega_lambda = this.dark_energy.parameter_value;
        let omega_r = this.calcul_omega_r();
        let a = omega_r;
        let b = 4 * omega_r + omega_m;
        let C = 5 * omega_r + 2 * omega_m - omega_lambda + 1;
        let p = -3 * Math.pow(b / a, 2) / 8 + C / a;
        var Ua = -2 * p / 6;
        var Ub = 5e10;
        let ex = 0.0000001;
        let Dev_Ua = this.derive_function_u(Ua);
        let Dev_Ub = this.derive_function_u(Ub);
        if (Dev_Ua * Dev_Ub < 0) {
            let extremum = Number(this.dichotomie(Ua, Ub, this.derive_function_u, 0, ex));
            if (this.function_u(extremum) > 0) {
                if (this.function_u(0) * this.function_u(Ub) < 0) {
                    let u_root = Number(this.dichotomie(0, Ub, this.function_u, 0, ex));
                    return u_root;
                }
                else if (this.function_u(0) * this.function_u(-Ub) < 0) {
                    let u_root = Number(this.dichotomie(-Ub, 0, this.function_u, 0, ex));
                    return u_root;
                }
                else {
                    return NaN;
                }
            }
            else if (this.function_u(extremum) == 0) {
                return extremum;
            }
            else {
                let u_root = Number(this.dichotomie(extremum, Ub, this.function_u, 0, ex));
                return u_root;
            }
        }
        else {
            if (this.function_u(0) * this.function_u(Ub) < 0) {
                let u_root = Number(this.dichotomie(0, Ub, this.function_u, 0, ex));
                return u_root;
            }
            else if (this.function_u(0) * this.function_u(-Ub) < 0) {
                let u_root = Number(this.dichotomie(-Ub, 0, this.function_u, 0, ex));
                return u_root;
            }
            else {
                return NaN;
            }
        }
    }
    function_u(u) {
        let omega_m = this._matter_parameter;
        let omega_lambda = this.dark_energy.parameter_value;
        let omega_r = this.calcul_omega_r();
        let a = omega_r;
        let b = 4 * omega_r + omega_m;
        let C = 5 * omega_r + 2 * omega_m - omega_lambda + 1;
        let d = 2 * omega_r + omega_m - 2 * omega_lambda + 2;
        let e = 1;
        let p = -3 * Math.pow(b / a, 2) / 8 + C / a;
        let q = Math.pow(b / a, 3) / 8 - (C / a) * (b / a) / 2 + d / a;
        let r = -3 * Math.pow(b / a, 4) / 256 + (C / a) * Math.pow(b / a, 2) / 16 - (d / a) * (b / a) / 4 + e / a;
        return (Math.pow(u, 3) - p * Math.pow(u, 2) - 4 * r * u + (4 * p * r - q * q));
    }
    derive_function_u(u) {
        let omega_m = this._matter_parameter;
        let omega_lambda = this.dark_energy.parameter_value;
        let omega_r = this.calcul_omega_r();
        let a = omega_r;
        let b = 4 * omega_r + omega_m;
        let C = 5 * omega_r + 2 * omega_m - omega_lambda + 1;
        let d = 2 * omega_r + omega_m - 2 * omega_lambda + 2;
        let e = 1;
        let p = -3 * Math.pow(b / a, 2) / 8 + C / a;
        let r = -3 * Math.pow(b / a, 4) / 256 + (C / a) * Math.pow(b / a, 2) / 16 - (d / a) * (b / a) / 4 + e / a;
        return (3 * u * u - 2 * p * u - 4 * r);
    }
    dichotomie(BornInf, BornSup, fonction, cible, ex) {
        let z_inf = BornInf;
        let z_sup = BornSup;
        let eps = ex;
        let dm_z_inf = fonction(this, z_inf);
        let dm_z_sup = fonction(this, z_sup);
        let max_iterations = 500;
        let j = 0;
        while (j < 500) {
            let zc = (z_inf + z_sup) / 2.0;
            let dm_zc = fonction(this, zc);
            if (((z_sup - z_inf) / 2) < ex) {
                if (Math.abs(zc / ex) < 100) {
                    ex = ex * 1e-5;
                }
                else {
                    return zc.toExponential(3);
                }
            }
            else if (isNaN(dm_zc)) {
                return NaN;
            }
            else if ((dm_zc - cible) * (dm_z_sup - cible) < 0) {
                z_inf = zc;
                dm_z_inf = dm_zc;
                j += 1;
            }
            else {
                z_sup = zc;
                dm_z_sup = dm_zc;
                j += 1;
            }
        }
    }
    metric_distance_simpson(Simu, z) {
        let distance;
        let curvature = Simu.calcul_omega_k();
        distance = Simu.simpson(Simu, Simu.integral_distance, 0, Number(z), 100000);
        if (curvature < 0) {
            distance =
                Math.sinh(Math.sqrt(Math.abs(curvature)) * distance) /
                    Math.sqrt(Math.abs(curvature));
        }
        else if (curvature > 0) {
            distance =
                Math.sin(Math.sqrt(Math.abs(curvature)) * distance) /
                    Math.sqrt(Math.abs(curvature));
        }
        distance *= Simu.constants.c / Simu._H0parsec;
        return distance;
    }
    /**
     * Will be used for the reverse calculation of the shift
     */
    Integral_dm(x) {
        return Math.sqrt(Math.abs(this.calcul_omega_k())) *
            this.simpson(this, this.integral_distance, 0, Number(x), 100000);
    }
}

let s = new Simulation_universe();


let t0=25;
let h0=50;
let omega_m=0.8;
let omega_lambda=0.1;
s.hubble_cst=h0;
s.temperature=t0;
s.matter_parameter=omega_m;
s.dark_energy.parameter_value=omega_lambda;



console.log(s.calcul_omega_k()+s.calcul_omega_r()+s.matter_parameter
+s.dark_energy.parameter_value);
s.single_fluid("radiation");

console.log("omega_k : ",s.calcul_omega_k());
console.log("omega_r : ",s.calcul_omega_r());
console.log("omega_m : ",s.matter_parameter);
console.log("omega_lamda : ",s.dark_energy.parameter_value);
console.log("age univ : ",s.universe_age());


console.log("H0 : ",s.H0parsec);

console.log("c : ",s.constants.G);