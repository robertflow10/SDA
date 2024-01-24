
final public  class aerosol {

    static double lm=2.0, gm=0.3, k=0.00001;

    static public         void main(String[] args) {


        int i, j;
        int a;     // user defined (argv[1]) total number of gas molecules in simulation
        int k1, ts; // for time stepping, including user defined (argv[2]) number of timesteps to integrate
        int rc;      // return code

        if       (args.length > 0)
        {         a = Integer.parseInt(args[0]);
            System.out.println("Number of gas molecules in the simulation: "  + Double.toString(a));
        }
        else a = 100;
        if      (args.length > 1)
        {
            ts = Integer.parseInt(args[1]);
            System.out.println("Number of time steps in the simulation: "  + Double.toString(ts));
        } else ts =  25;

        System.out.printf("Initializing for %d particles in x,y,z space...%n", a);

        // 1D array for mass
        double[] M
                = new double[a];
        // 1D array for (x,y,z) position
        double[] x
                = new double[a];
        double[] y
                = new double[a];
        double[] z
                = new double[a];
        //1D array for (vx, vy, vz) velocity
        double[] vx
                = new double[a];
        double[] vy
                = new double[a];
        double[] vz
                = new double[a];

        double dx, dy, dz, d, F, GRAVCONST=0.001, b=300;
        double ax, ay, az;

        // 1D array for each particle's component that will evaporate
        double[] gas = new double[a];
        double[] liquid = new double[a];
        double[] loss_rate = new double[a];

        // save previous values whilst doing global updates
        double[] o1 = new double[a];
        double[] o2 = new double[a];
        double[] o3 = new double[a];
        double[] om = new double[a];
        double totalmass, E;  // for stats

        double t1 = System.currentTimeMillis();

        // initialise
        rc = INIT(M, x, y, z, vx, vy, vz, gas, liquid, loss_rate, a);
        if (rc != 0) {
            System.out.println("Error in INIT");
            System.exit(1);
        }
        else System.out.println("INIT COMPLETE");

        totalmass = 0.0;
        for (i=0; i<a; i++)


        {M[i] = gas[i]*gm + liquid[i]*lm;
            totalmass += M[i];        }

        E = calc_system_energy(totalmass, vx, vy, vz, a);

        System.out.printf("Time 0. System energy=%g%n", E);


  /*
     MAIN TIME STEPPING LOOP

     We 'save' old (entry to timestep loop) values to use on RHS of:

     For each aerosol particle we will: calc forces due to other
     particles & update change in velocity and thus position; then we
     condense some of the gas to liquid which changes the mass of the
     particle so we then determine its new velocity via conservation
     of kinetic energy.

     Having looped over the particles (using 'old' values on the right
     hand side as a crude approximation to real life) we then have
     updated all particles independently.

     We then determine the total mass & the system energy of the
     system.

     The final statement of each time-stepping loop is to decrease the
     temperature of the system.

     After completing all time-steps, we output the the time taken and
     the centre of mass of the final system configuration.

  */

        System.out.printf("Now to integrate for %d timesteps%n", ts);

        // time=0 was initial conditions
        for (k1=1; k1<=ts; k1++) {

            // LOOP1: take snapshot to use on RHS when looping for updates
            for (i=0; i<a; i++) {
                o1[i] = x[i];
                o2[i] = y[i];
                o3[i] = z[i];
                om[i] = M[i];
            }

            double temp_d, temp_z;

            // LOOP2: update position etc per particle (based on old data)
            for(i=0; i<a; i++) {
                // calc forces on body i due to particles (j != i)
                for (j=0; j<a; j++) {
                    if (j != i) {
                        dx = o1[j] - x[i];
                        dy = o2[j] - y[i];
                        dz = o3[j] - z[i];
                        d = Math.max(Math.sqrt(dx * dx + dy * dy + dz * dz), 0.01);
                        F = GRAVCONST * M[i] * om[j] / (d*d);
                        // calculate acceleration due to the force, F
                        ax = (F/M[i]) * dx/d;
                        ay = (F/M[i]) * dy/d;
                        az = (F/M[i]) * dz/d;
                        // approximate velocities in "unit time"
                        vx[i] += ax;
                        vy[i] += ay;
                        vz[i] += az;
                    }
                }
                // calc new position
                x[i] = o1[i] + vx[i];
                y[i] = o2[i] + vy[i];
                z[i] = o3[i] + vz[i];
                // temp-dependent condensation from gas to liquid
                gas[i] *= loss_rate[i] * Math.exp(-k*b);
                liquid[i] = 1.0 - gas[i];
                M[i] = gas[i]*gm + liquid[i]*lm;
                // conserve energy systems 0.5*m*v*v remains constant
                double v_squared = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
                double factor = Math.sqrt(om[i]*v_squared/M[i])/Math.sqrt(v_squared);
                vx[i] *= factor;
                vy[i] *= factor;
                vz[i] *= factor;
            } // end of LOOP 2
            _a(x,y,z, vx,vy,vz, gas, liquid, a);
            totalmass = 0.0;
            for (i=0; i<a; i++) {
                totalmass += M[i];
            }
            E = calc_system_energy(totalmass, vx, vy, vz, a);
            System.out.printf("At end of timestep %d with temp %f the system energy=%g and total aerosol mass=%g%n", k1, b, E, totalmass);
            // temperature drops per timestep
            b *= 0.99999;
        } // time steps

        System.out.printf("Time to init+solve %d molecules for %d timesteps is %g ms%n", a, ts, System.currentTimeMillis()-t1);
        // output a metric (centre of mass) for checking
        double[]com  = new double[3];
        calcCentre_Mass(com, x,y,z,M,totalmass,a);
        System.out.printf("Centre of mass = (%g,%g,%g)%n", com[0], com[1], com[2]);
        System.exit(0);
    } // main

    private static void _a(double[] x, double[] y, double[] z, double[] vx, double[] vy, double[] vz, double[] gas, double[] liquid, int n) {
        int i;
        System.out.printf("num \t position (x,y,z) \t velocity (vx, vy, vz) \t mass (gas, liquid)%n");
        for (i=0; i<n; i++) {
            System.out.printf("%d \t %f %f %f \t %f %f %f \t %f %f\n", i, x[i], y[i], z[i], vx[i], vy[i], vz[i], gas[i]*gm, liquid[i]*lm);
        }
    }


    private static double calc_system_energy(double M, double[] vx, double[] vy, double[] vz, int n) {
        double TOTAL_ENERGY = 0.0;
        for (int i = 0; i < n; i++) {
            TOTAL_ENERGY += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
        }
        TOTAL_ENERGY = 0.5 * M * TOTAL_ENERGY;
        return TOTAL_ENERGY / (double) n;
    }


    private static void calcCentre_Mass(double[] com, double[] x, double[] y, double[] z, double[] mass, double totalMass, int N) {
        int i, axis;
        // calculate the centre of mass, com(x,y,z)
        for (axis=0; axis<2; axis++) {
            com[0] = 0.0;     com[1] = 0.0;     com[2] = 0.0;
            for (i=0; i<N; i++) {
                com[0] += mass[i]*x[i];
                com[1] += mass[i]*y[i];
                com[2] += mass[i]*z[i];
            }
            com[0] /= totalMass;
            com[1] /= totalMass;
            com[2] /= totalMass;
        }
    }


    private static int INIT(double[] mass, double[] x, double[] y, double[] z, double[] vx, double[] vy, double[] vz, double[] gas, double[] liquid, double[] loss_rate, int a) {
        // random numbers to set initial conditions - do not parallelise or amend order of random number usage
        int i;
        double comp;
        double min_pos = -50.0, mult = +100.0, maxVel = +10.0;

        for (i=0; i<a; i++) {
            x[i] = min_pos + mult*Math.random();
            y[i] = min_pos + mult*Math.random();
            z[i] = 0.0 + mult*Math.random();
            vx[i] = -maxVel + 2.0*maxVel*mult*Math.random();
            vy[i] = -maxVel + 2.0*maxVel*mult*Math.random();
            vz[i] = -maxVel + 2.0*maxVel*mult*Math.random();
            // proportion of aerosol that evaporates
            comp = .5 + Math.random()/2.0;
            loss_rate[i] = 1.0 - Math.random()/25.0;
            // aerosol is component of gas and (1-comp) of liquid
            gas[i] = comp;
            liquid[i] = (1.0-comp);
        }
        return 0;
    } // init


}
