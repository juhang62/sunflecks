/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package sunflecks;

import static java.lang.Math.*;

//comment which class from which package
/**
 *
 * @author lfhan
 */
public final class Sunflecks {

    //stomata para
    public double gsmax = 1.8000e-01;
    public double gsmin = 9.1777e-02;
    public double gb = 1440000;
    public double alphag = 2.0866e-03;
    public double thetag = 9.0867e-01;
    /**/
    public double taugi = 2.0867e+02;
    public double taugd = 3.6844e+02;
    public double taupi = 1.8467e+02;
    public double tauw = 1.8519e+02;
    /**/
    public double rd = 1.7; //miumol
    public double vjmax = 30.0; //miumol  = 2 * this.vcmax;
    public double vcmax = 98.0; //miumol
    public double vfmax = 2.4*98; //this.vfmax = 1.8 * this.vcmax;
    public double vfmin = 0.1;
    public double vcmin = 0.01;
    /**/
    public double alphac = 6.0088e-03; //est.
    public double thetac = 3.5740e-01;
    public double alphaf = 6.0088e-03;
    public double thetaf = 3.5740e-01;
    public double alphaj = 6.0088e-03;
    public double thetaj = 3.5740e-01;
    /**/
    public double rmax = 196; //miu mol, also used in dy function
    public double tmax = 196; //miu mol, also used in dy function
    public double kr = 2; //need est. miu mol, also used in dy function
    public double kt = 2;  //need est. miu mol, also used in dy function 
    public double gamma = 4.44;//also used in dy function
    public double ka = 1;
    public double kc = 31;
    public double ko = 15500;
    public double patm = 100000;
    public double ca = 39;//changed
    public double po2 = 21000;
    //psyn timeconst
    public double psi = 0.04;
    public double tauci = 50;
    public double taucd = 120;
    public double taufi = 30;
    public double taufd = 50;
    //State variables
    public double signal;
    public double pi;
    public double water;
    public double poolT;
    public double poolR;
    public double poolG;
    public double ci;
    public double irr;
    public double gs;
    public double vj;
    public double vf;
    public double vc;
    public double ass;
    public double seq;
    public double curtime;
    public double dt = 0.1;
    public boolean laststep;
    public boolean ifInputgs = false;
    private boolean ifloop = false;

    public static void main(String[] args) {
        Sunflecks obj = new Sunflecks(50);
        double[] time = new double[3600];
        double[] irrin = new double[3600];
        double[] gsco2 = new double[3600];
        double[] ambco2 = new double[3600];
        for (int i = 0; i < 2; i++) {
            time[i] = i + 1;
            irrin[i] = 50;
//            gsco2[i] = 0.03;
//            ambco2[i] = 39;

        }
        for (int i = 2; i < 3600; i++) {
            time[i] = i + 1;
            irrin[i] = 1300;
//            gsco2[i] = 0.03;
//            ambco2[i] = 39;
        }
//        double[] time={1,2,3,4,5};
//        double[] irr={100,200,300,400,500};
        //double[] ass = obj.caliDynA(time, irr, gsco2, ambco2);
        double[][] result = obj.rawrun(time, irrin);
        System.out.println(result[0][1]);
        for (int i = 0; i < 5; i++) {
        System.out.println(result[0][i]);
        }    
    }

    public double[][] rawrun(double[] time, double[] irr) {
        //double dt=0.01;
        int nDatapts = time.length;
        double timeintv;
        double[] assco2 = new double[nDatapts];
        double[] vcout = new double[nDatapts];
        double[] vfout = new double[nDatapts];
        double[] vjout = new double[nDatapts];
        double[] gsout = new double[nDatapts];
        //Sunflecks simSunflecks=new Sunflecks(timeirr[0][1]); 
        this.calcDyn(irr[0], time[0]);
        assco2[0] = this.ass;
        vcout[0] = this.vc;
        vfout[0] = this.vf;
        vjout[0] = this.vj;
        gsout[0] = this.gs;
        for (int i = 1; i < nDatapts; i++) {
            timeintv = time[i] - time[i - 1];
            this.calcDyn(irr[i], timeintv);
            //get An at each data point and write into output file
            assco2[i] = this.ass;
            vcout[i] = this.vc;
            vfout[i] = this.vf;
            vjout[i] = this.vj;
            gsout[i] = this.gs;
        }
        double[][] out = {assco2, vcout, vfout, vjout, gsout};
        return out;
    }

    public double[] caliDynA(double[] time, double[] irr, double[] gsco2, double[] ambco2) {
        int nDatapts = time.length;
        double timeintv;
        double[] assout = new double[nDatapts];
        this.gs = gsco2[0];
        this.ca = ambco2[0];
        this.calcSteady(irr[0]);
        this.calcDyn(irr[0], time[0]);
        assout[0] = this.ass;
        for (int i = 1; i < nDatapts; i++) {
            this.gs = gsco2[i];
            this.ca = ambco2[i];
            timeintv = time[i] - time[i - 1];
            this.calcDyn(irr[i], timeintv);
            assout[i] = this.ass;
        }
        return assout;
    }

    //constructor be called when instantiate the class
    public Sunflecks(double iniIRR) {
        this.calcSteady(iniIRR);
    }

    private String[] getOutput(Double curtime) {
        String[] strArr = new String[15];
        strArr[0] = String.valueOf(curtime);
        strArr[1] = String.valueOf(this.ass);
        strArr[2] = String.valueOf(this.irr);
        strArr[3] = String.valueOf(this.ci);
        strArr[4] = String.valueOf(this.vc);
        strArr[5] = String.valueOf(this.vf);
        strArr[6] = String.valueOf(this.vj);
        strArr[7] = String.valueOf(this.poolT);
        strArr[8] = String.valueOf(this.poolR);
        strArr[9] = String.valueOf(this.poolG);
        strArr[10] = String.valueOf(this.gs);
        strArr[11] = String.valueOf(this.signal);
        strArr[12] = String.valueOf(this.pi);
        strArr[13] = String.valueOf(this.water);
        strArr[14] = String.valueOf(this.seq);
        return strArr;
    }

    //private ici does not converge 315, 1242..
    public void calcSteady(double irr) {
        this.irr = irr;
        this.vj = getvxeq(this.alphaj, this.thetaj, this.vjmax, 0);
        this.vf = getvxeq(this.alphaf, this.thetaf, this.vfmax, this.vfmin);
        vc = getvxeq(this.alphac, this.thetac, this.vcmax, this.vcmin);
        if (ifInputgs == false) {
            updateStdygs();
        }
        double newci = 38.9; //first guess
        this.ci = 39;
        int ici = 0;
        while (abs(newci - this.ci) > 0.0001 && ici < 1000) {
            ici = ici + 1;
//            //double oldci=this.ci;
//            this.ci= newci;
//            double fx1=fci(this.ci);
//            double fx2=fci(this.ci-0.01);
//            //Netwon's Method: xnew=xold-f(xold)/f'(xold)
//            newci=this.ci-fx1*(0.01/(fx1-fx2)); 
            double oldci = this.ci;
            this.ci = newci;
            double fx1 = fci(this.ci);
            double fx2 = fci(oldci);
            //secent method
            newci = this.ci - fx1 * ((this.ci - oldci) / (fx1 - fx2));
        }
        if (ici == 999) {
            System.out.println("ici error");
            System.out.println(irr);
        }
        double gg = getgg();
        this.ass = (this.ca - this.ci) * gg / this.patm;
    }

    public void calcDyn(double irr, double timeinv) {
        this.irr = irr;
        this.curtime = 0.0;
        //this.laststep = false;
        while (this.curtime < timeinv) {
            if (this.curtime + dt > timeinv) {
                this.dt = timeinv - this.curtime;
//                this.laststep = true;
//            } else {
//                this.laststep = false;
            }
            //update biochem first in order fix gs bug?
            this.vj = getvxeq(this.alphaj, this.thetaj, this.vjmax, 0);
            double vfeq = getvxeq(this.alphaf, this.thetaf, this.vfmax, this.vfmin);
            double vceq = getvxeq(this.alphac, this.thetac, this.vcmax, this.vcmin);
            double wc = getwc(this.ci);
            double[] tau = {this.tauci, this.taucd, this.taufi, this.taufd}; //time const
            //create a object of Differential Equation using my class Dfun
            Dfun photoDif = new Dfun(tau, this.rmax, this.tmax, this.kr,
                    this.kt, this.gamma, this.psi, this.vj, vfeq, vceq, this.ci, wc);
            //initial value of state variables
            double[] y0 = {this.vc, this.vf, this.poolT, this.poolR, this.poolG};
            //calculate those values after one time step
            //double[] yf=rk4(y0,dt,photoDif);
            double[] yf = rkf45beta(y0, photoDif); //use my rkf45 method 
            //important to sync with gs time step
            if (ifloop == true) {
                continue;
            }
            this.vc = yf[0];
            this.vf = yf[1];
            this.poolT = yf[2];
            this.poolR = yf[3];
            this.poolG = yf[4];
            wc = getwc(this.ci);
            this.ass = getass1(wc);

            //update gs, solve for ci
            if (ifInputgs == false) {
                updateDyngs(dt);
            }
            double gg = getgg();
            this.ci = this.ca - ((this.ass * this.patm) / gg);
            double newci = this.ci + 1; //first guess
            for (int ci_i = 1; ci_i < 500; ci_i++) {
                double oldci = this.ci;
                this.ci = newci;
                double fx1 = dynfci(this.ci);
                double fx2 = dynfci(oldci);
                //secent method
                newci = this.ci - fx1 * ((this.ci - oldci) / (fx1 - fx2));
                if (abs(newci - this.ci) < 0.0001) {
                    break;
                }
            }
            if (abs(newci - this.ci) > 0.0001) {
                System.out.println("dyn ci error");
            }
            
            wc = getwc(this.ci);
            this.ass = getass1(wc);
        }
        //System.out.println(this.ass);
    }

    //dynamic function of ci which we want it to be zero
    public double dynfci(double localci) {
        double wc = getwc(localci);
        double ass1 = getass1(wc);
        double gg = getgg();
        double ass2 = (this.ca - localci) * gg / this.patm;
        return ass1 - ass2;
    }

    //steady state function of ci which we want it to be zero
    public double fci(double localci) {
        double wc = getwc(localci);
        this.poolT = this.tmax;
        this.poolR = getR();
        double eqn2 = geteqn2(localci, wc);    //eqn2 in Pearcy1997 
        //bisection method 
        //first find T1 and T2 which statify eqn2(T1)<0 , eqn2(T2)>0 
        double T2 = 0.0;
        double T1 = 0.0;
        if (eqn2 > 0) {
            T2 = this.tmax;
            //find t1
            for (int findt1_i = 50; findt1_i < 5000; findt1_i += 100) {
                double Tportion = T2 / (findt1_i + 1);
                for (int findt1_j = 1; findt1_j <= findt1_i; findt1_j++) {
                    T1 = T2 - Tportion * findt1_j;
                    this.poolT = T1;
                    this.poolR = getR();
                    eqn2 = geteqn2(localci, wc);
                    if (eqn2 < 0) {
                        break;
                    }
                }
                if (eqn2 < 0) {
                    break;
                }
            }
            if (eqn2 > 0) {
                System.out.println("eqn2 error");
            }
        } else { //eqn2<0
            T1 = this.tmax;
            //find t2
            for (int findt2_i = 50; findt2_i < 5000; findt2_i += 100) {
                double Tportion = T1 / (findt2_i + 1);
                for (int findt2_j = 1; findt2_j <= findt2_i; findt2_j++) {
                    T2 = T1 - Tportion * findt2_j;
                    this.poolT = T2;
                    this.poolR = getR();
                    eqn2 = geteqn2(localci, wc);
                    if (eqn2 > 0) {
                        break;
                    }
                }
                if (eqn2 > 0) {
                    break;
                }
            }
            if (eqn2 < 0) {
                System.out.println("eqn2 error");
            }
        }

        //bisection
        for (int bisec_i = 1; bisec_i < 1000; bisec_i++) {
            if (abs(eqn2) < 0.001) {
                break;
            }
            this.poolT = (T1 + T2) / 2.0;
            this.poolR = getR();
            eqn2 = geteqn2(localci, wc);
            if (eqn2 > 0) {
                T2 = this.poolT;
            } else {
                T1 = this.poolT;
            }
        }
        if (abs(eqn2) > 0.001) {
            System.out.println("bisec error");
        }
        this.poolG = ((2.0 * this.gamma / localci) * wc
                * (this.poolR / (this.kr + this.poolR))) / this.psi;
        double ass1 = getass1(wc);
        double gg = getgg();
        double ass2 = (this.ca - localci) * gg / this.patm;
        return ass1 - ass2;
    }

    //method to update state varibles of steady state stomata
    public void updateStdygs() {
        this.seq = getSeq();
        this.gs = this.seq * this.gsmax;
        this.signal = this.seq;
        this.pi = this.seq;
        this.water = this.seq;
    }

    //method to update state varibles of dynamic stomata
    public void updateDyngs(double dt) {
        this.seq = getSeq();
        Dfun gsDif = new Dfun(this.taugi, this.taugd, this.taupi, this.tauw, this.seq);
        double[] y0 = {this.signal, this.pi, this.water};
        double[] yf = rk4(y0, dt, gsDif);
        this.signal = yf[0];
        this.pi = yf[1];
        this.water = yf[2];
        this.gs = this.gsmax * this.water;
    }

    //runge-kutta
    public double[] rk4(double[] y0, double dt, Dfun dfun) {
        double[] k1 = dfun.getdf(y0);
        k1 = scaleArray(k1, dt);
        double[] y1 = addArray(y0, scaleArray(k1, 0.5));
        double[] k2 = dfun.getdf(y1);
        k2 = scaleArray(k2, dt);
        double[] y2 = addArray(y0, scaleArray(k2, 0.5));
        double[] k3 = dfun.getdf(y2);
        k3 = scaleArray(k3, dt);
        double[] y3 = addArray(y0, k3);
        double[] k4 = dfun.getdf(y3);
        k4 = scaleArray(k4, dt);
        double[] sumk = addMultiArray(k1, scaleArray(k2, 2), scaleArray(k3, 2), k4);
        sumk = scaleArray(sumk, 1.0 / 6);
        double[] yf = addArray(y0, sumk);
        return yf;
    }

    public double[] rkf45(double[] y0, Dfun dfun, double tspan) {
        double dt = 0.01;
        double dtmax = 0.5;
        double dtmin = 0.0005;
        double t = 0;
        double reltol = 1e-3;
        double abstol = 1e-3;
        double tol = 1;
        int nvar = y0.length;
        double[][] rkfinmd = new double[5][5];
        rkfinmd[0][0] = 1.0 / 4;
        rkfinmd[1][0] = 3.0 / 32;
        rkfinmd[1][1] = 9.0 / 32;
        rkfinmd[2][0] = 1932.0 / 2197;
        rkfinmd[2][1] = -7200.0 / 2197;
        rkfinmd[2][3] = 7296.0 / 2197;
        rkfinmd[3][0] = 439.0 / 216;
        rkfinmd[3][1] = -8.0;
        rkfinmd[3][2] = 3680.0 / 513;
        rkfinmd[3][3] = -845.0 / 4104;
        rkfinmd[4][0] = -8.0 / 27;
        rkfinmd[4][1] = 2.0;
        rkfinmd[4][2] = -3544.0 / 2565;
        rkfinmd[4][3] = 1859.0 / 4104;
        rkfinmd[4][4] = -11.0 / 40;

        double[][] rkffi = new double[2][6];
        rkffi[0][0] = 25.0 / 216;
        rkffi[0][1] = 0.0;
        rkffi[0][2] = 1408.0 / 2565;
        rkffi[0][3] = 2197.0 / 4104;
        rkffi[0][4] = -1.0 / 5;
        rkffi[0][5] = 0.0;
        rkffi[1][0] = 16.0 / 135;
        rkffi[1][1] = 0.0;
        rkffi[1][2] = 6656.0 / 12825;
        rkffi[1][3] = 28561.0 / 56430;
        rkffi[1][4] = -9.0 / 50;
        rkffi[1][5] = 2.0 / 55;
        boolean laststep;
        while (t < tspan) {
            if (t + dt > tspan) {
                dt = tspan - t;
                laststep = true;
            } else {
                laststep = false;
            }
            double[][] k = new double[6][nvar];
            k[0] = dfun.getdf(y0);
            k[0] = scaleArray(k[1], dt);
            double[] y = addArray(y0, scaleArray(k[0], rkfinmd[0][0]));
            for (int i = 1; i < 6; i++) {
                k[i] = dfun.getdf(y);
                k[i] = scaleArray(k[i], dt);
                if (i == 5) {
                    break;
                }
                y = addArray(y0, scaleArray(k[0], rkfinmd[i][0]));
                for (int j = 1; j <= i; j++) {
                    y = addArray(y, scaleArray(k[j], rkfinmd[i][j]));
                }
            }

            //double[] yf4 = new double[nvar]; //caution:pass by value; no yf4=y0
            //double[] yf5 = new double[nvar];
            double[] yf4 = addArray(y0, scaleArray(k[0], rkffi[0][0]));
            double[] yf5 = addArray(y0, scaleArray(k[0], rkffi[1][0]));
            for (int i = 1; i < 6; i++) {
                yf4 = addArray(yf4, scaleArray(k[i], rkffi[0][i]));
                yf5 = addArray(yf5, scaleArray(k[i], rkffi[1][i]));
            }
            double err = normmax(addArray(yf4, scaleArray(yf5, -1)));
            double ratio = err / (reltol * minabs(yf4) + abstol);
            //double ratio=err/dt;
            //double theta=0.84*pow((reltol*normmax(yf4)/ratio),0.25);
            if (ratio < 1) {
                t = t + dt;
                y0 = yf4;
                //dt=theta*dt;
                dt = 0.9 * dt * pow(ratio, -0.2);
                //dt=dt*0.840896*pow((dt/err),0.25);
            } else {
                //dt=theta*dt;
                dt = 0.9 * dt * pow(ratio, -0.25);
                //dt=dt*0.840896*pow((dt/err),0.25);
            }

            if (dtmax < dt) {
                dt = dtmax;
            }
            if (dtmin > dt) {
                dt = 0.001;
                if (laststep == false) {
                    System.out.println("stiff");
                }
            }
        }

        return y0;
    }

    public double[] rkf45beta(double[] y0, Dfun dfun) {
        double reltol = 1e-3;
        double abstol = 1e-3;
        int nvar = y0.length;
        double[][] rkfinmd = new double[5][5];
        rkfinmd[0][0] = 1.0 / 4;
        rkfinmd[1][0] = 3.0 / 32;
        rkfinmd[1][1] = 9.0 / 32;
        rkfinmd[2][0] = 1932.0 / 2197;
        rkfinmd[2][1] = -7200.0 / 2197;
        rkfinmd[2][3] = 7296.0 / 2197;
        rkfinmd[3][0] = 439.0 / 216;
        rkfinmd[3][1] = -8.0;
        rkfinmd[3][2] = 3680.0 / 513;
        rkfinmd[3][3] = -845.0 / 4104;
        rkfinmd[4][0] = -8.0 / 27;
        rkfinmd[4][1] = 2.0;
        rkfinmd[4][2] = -3544.0 / 2565;
        rkfinmd[4][3] = 1859.0 / 4104;
        rkfinmd[4][4] = -11.0 / 40;

        double[][] rkffi = new double[2][6];
        rkffi[0][0] = 25.0 / 216;
        rkffi[0][1] = 0.0;
        rkffi[0][2] = 1408.0 / 2565;
        rkffi[0][3] = 2197.0 / 4104;
        rkffi[0][4] = -1.0 / 5;
        rkffi[0][5] = 0.0;
        rkffi[1][0] = 16.0 / 135;
        rkffi[1][1] = 0.0;
        rkffi[1][2] = 6656.0 / 12825;
        rkffi[1][3] = 28561.0 / 56430;
        rkffi[1][4] = -9.0 / 50;
        rkffi[1][5] = 2.0 / 55;

        double[][] k = new double[6][nvar];
        k[0] = dfun.getdf(y0);
        k[0] = scaleArray(k[1], dt);
        double[] y = addArray(y0, scaleArray(k[0], rkfinmd[0][0]));
        for (int i = 1; i < 6; i++) {
            k[i] = dfun.getdf(y);
            k[i] = scaleArray(k[i], dt);
            if (i == 5) {
                break;
            }
            y = addArray(y0, scaleArray(k[0], rkfinmd[i][0]));
            for (int j = 1; j <= i; j++) {
                y = addArray(y, scaleArray(k[j], rkfinmd[i][j]));
            }
        }

        //double[] yf4 = new double[nvar]; //caution:pass by value; no yf4=y0
        //double[] yf5 = new double[nvar];
        double[] yf4 = addArray(y0, scaleArray(k[0], rkffi[0][0]));
        double[] yf5 = addArray(y0, scaleArray(k[0], rkffi[1][0]));
        for (int i = 1; i < 6; i++) {
            yf4 = addArray(yf4, scaleArray(k[i], rkffi[0][i]));
            yf5 = addArray(yf5, scaleArray(k[i], rkffi[1][i]));
        }
        double err = normmax(addArray(yf4, scaleArray(yf5, -1)));
        double ratio = err / (reltol * minabs(yf4) + abstol);
        //double ratio=err/dt;
        //double theta=0.84*pow((reltol*normmax(yf4)/ratio),0.25);
        if (ratio < 0.003) {
            this.curtime = this.curtime + this.dt;
            y0 = yf4;
            this.dt = 0.1;
            ifloop = false;
        } else if (ratio < 1) {
            this.curtime = this.curtime + this.dt;
            y0 = yf4;
            this.dt = 0.01;
            ifloop = false;
        } else if (dt <= 0.001) {
            this.curtime = this.curtime + this.dt;
            y0 = yf4;
            ifloop = false;
        } else {
            //dt=theta*dt;
            ifloop = true;
            this.dt = 0.001;
            //dt=dt*0.840896*pow((dt/err),0.25);
        }
        return y0;
    }

    public static double normmax(double[] vec) {
        double max = 0;
        for (int i = 0; i < vec.length; i++) {
            if (abs(vec[i]) > max) {
                max = abs(vec[i]);
            }
        }
        return max;
    }

    public static double minabs(double[] vec) {
        double min = 10000;
        for (int i = 0; i < vec.length; i++) {
            if (abs(vec[i]) < min) {
                min = abs(vec[i]);
            }
        }
        return min;
    }

    public static double[] scaleArray(double[] arr, double scaleFactor) {
        double[] result;
        result = new double[arr.length];
        for (int i = 0; i < arr.length; i++) {
            result[i] = arr[i] * scaleFactor;
        }
        return result;
    }

    public static double[] addArray(double[] arr1, double[] arr2) {
        int nelem = arr1.length;
        double[] result = new double[nelem];
        for (int i = 0; i < arr1.length; i++) {
            result[i] = arr1[i] + arr2[i];
        }
        return result;
    }

    public static double[] addMultiArray(double[]... arr) {
        int narr = arr.length;
        int nelem = arr[0].length;
        double[] sum = new double[nelem];
        for (int i = 0; i < narr; i++) {
            sum = addArray(sum, arr[i]);
        }
        return sum;
    }

    //get equilibrium signal 
    public double getSeq() {
        //double smin = this.gsmin / this.gsmax;
        double seqval = ((1 + this.gsmin + this.alphag * this.irr)
                - sqrt(pow(1 + this.gsmin + this.alphag * this.irr, 2)
                - 4 * this.thetag * (this.gsmin + this.alphag * this.irr))) / (2.0 * this.thetag);
        return seqval;
    }

    //get equilibrium vf, vj, vc
    private double getvxeq(double alphax, double thetax, double vmax, double vxmin) {
        double fn = (((alphax * this.irr + 1.0 - vxmin)
                - sqrt(pow(alphax * this.irr + 1.0 - vxmin, 2)
                - 4.0 * alphax * this.irr * thetax * (1.0 - vxmin)))
                / (2.0 * thetax)) + vxmin;
        double vxeq = vmax * fn;
        return vxeq;
    }

    private double getR() {
        double shit1 = (this.vj * (1 - this.poolT / this.tmax))
                * ((this.kt + this.poolT) / (this.vf * this.poolT));
        //no 5/3 in units of RuBP, corrected
        double pooR = (1 - shit1) * this.rmax;
        return pooR;
    }

    private double geteqn2(double localci, double wc) {
        double eqn2 = (((this.vf * this.poolT) / (this.kt + this.poolT))
                * (1 - this.poolR / this.rmax))
                - (1 + 2.0 * this.gamma / localci) * wc
                * (this.poolR / (this.kr + this.poolR));
        return eqn2;
    }

    private double getgg() {
        //convert units
        double gg = (this.gs * 1000000 * this.gb) / (this.gs * 1000000 + this.gb); 
        return gg;
    }

    private double getwc(double localci) {
        double wc = (this.vc * pow(localci, 2.0)) / ((localci + this.ka)
                * (localci + this.kc * (1 + this.po2 / this.ko)));
        return wc;
    }

    private double getass1(double wc) {
        double ass1 = wc * (this.poolR / (this.poolR + this.kr))
                - 0.5 * this.psi * this.poolG - this.rd;
        return ass1;
    }
}
