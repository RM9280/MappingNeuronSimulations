package src;

//Updated by Rebecca Mantione, orginal code from Chris Fietkiewicz. Hodgkin-Huxley model of a neuron.
// NOTE: Requires SwingGraphics.java for graphing.
import java.awt.*;
import java.util.ArrayList;
import java.util.Scanner;

public class HH_NeuronStarter3 {
	double dt;
// Constants
	double gBarK = 36; // mS/cm^2
	double gBarNa = 120; // mS/cm^2
	double gM = 0.3; // mS/cm^2
	double eK = -77; // mV
	double eNa = 50; // mV
	double vRest = -54.4; // mV
// Initial conditions
	double v = -65;
	double n = 0.3177;
	double m = 0.0529;
	double h = 0.5961;
	static ArrayList<Double> list;

// Constructor
	public HH_NeuronStarter3(double dt) {
		this.dt = dt;
	}

// Returns new membrane potential after a simulated time step.
// Receives a stimulus as a voltage which is added directly to membrane
//potential.
	public double calculateNextTimeStep(double stimulus) {
// Potassium current
		double alphan = 0.01 * (v + 55) / (1 - Math.exp(-(v + 55) / 10));
		double betan = 0.125 * Math.exp(-(v + 65) / 80);
// Sodium current
		double alpham = 0.1 * (v + 40) / (1 - Math.exp(-(v + 40) / 10));
		double betam = 4 * Math.exp(-(v + 65) / 18);
		double alphah = 0.07 * Math.exp(-(v + 65) / 20);
		double betah = 1 / (1 + Math.exp(-(v + 35) / 10));
// Differential equations
		double dn = (alphan * (1 - n) - betan * n);
		double dm = (alpham * (1 - m) - betam * m);
		double dh = (alphah * (1 - h) - betah * h);
		double dv = (-gBarNa * Math.pow(m, 3) * h * (v - eNa) - gBarK * Math.pow(n, 4) * (v - eK) - gM * (v - vRest)
				+ stimulus);
// Euler integration updates
		n = n + dn * dt;
		m = m + dm * dt;
		h = h + dh * dt;
		v = v + dv * dt;
		return v;
	}

	public double getdt() {
		return dt;
	}

	public static double simulate(double dt, int N, double[] t, double[] Vm, double stimAmplitude,
			double stimDuration) {
		HH_NeuronStarter n = new HH_NeuronStarter(dt); // Create a neuron
		list = new ArrayList<Double>();
		double firstMsec = 0.5; // Time for 1st stimulus
		int first = (int) Math.round(firstMsec / n.getdt());
		int intStimDuration = (int) Math.round(stimDuration / n.getdt()) - 1;
		double Vold = -65.0;
		double t1 = -1.0;
		double t2 = -1.0;
		System.out.println("Simulation settings: amplitude = " + stimAmplitude + ", dt = " + dt + " msec...");
		for (int i = 0; i < N; i++) {
			t[i] = i * n.getdt();
			if (i >= first && i <= (first + intStimDuration))
				Vm[i] = n.calculateNextTimeStep(stimAmplitude);
			else // Normally, recalcuate membrane potential with no stimulus
				Vm[i] = n.calculateNextTimeStep(0.0);

			if (Vold < 0 && Vm[i] >= 0) {
				list.add(t[i]);
			}
			Vold = Vm[i];
		}
		if (list.size() >= 2) {
			return 1000 / (list.get(list.size() - 1) - list.get(list.size() - 2));
		} else {
			return 0;
		}
	}

// Main method that runs the simulation
	public static void main(String[] args) {
		double[] arrayfreq = new double[22];
		double[] time = new double[22];
//		Scanner reader = new Scanner(System.in); // Reading from System.in
//		System.out.println("Enter starting amplitude: ");
		double n = 0; // Scans the next token of the input as an int
		// Once finished
//		reader.close();
//		System.out.println("readin" + n);
		double dt = 0.01; // Time step duration
		double duration = 50.0; // Simulation duration
		double stimDuration = 50.0; // Stimulus duration
// Experiment
		double stimAmplitude = n;
		int N = (int) (duration / dt);
		int round = 0;
		SwingGraphics grapher = new SwingGraphics(); // Begin graphing
		SwingGraphics grapher2 = new SwingGraphics(); // Begin graphing
		for (double i = n; i < n + 22; i++) {
			time[round] = i;
			stimAmplitude = i;
			double[] t1 = new double[N];
			double[] Vm1 = new double[N];
			double freq = simulate(dt, N, t1, Vm1, stimAmplitude, stimDuration);
			arrayfreq[(int) i] = freq;
			if (round == 1) {
				grapher.graph(t1, Vm1, Color.PINK);
			}
			if (round == 2) {
				grapher.graph(t1, Vm1, Color.ORANGE);
			}
			if (round == 3) {
				grapher.graph(t1, Vm1, Color.CYAN);
			}
			if (round == 4) {
				grapher.graph(t1, Vm1, Color.MAGENTA);
			}
			if (round == 5) {
				grapher.graph(t1, Vm1, Color.BLACK);
			}
			if (round == 6) {
				grapher.graph(t1, Vm1, Color.YELLOW);
			}
			if (round == 7) {
				grapher.graph(t1, Vm1, Color.RED);
			}
			if (round == 8) {
				grapher.graph(t1, Vm1, Color.GRAY);
			}
			if (round == 9) {
				grapher.graph(t1, Vm1, Color.LIGHT_GRAY);
			}
			if (round == 10) {
				grapher.graph(t1, Vm1, Color.black);
			}
			if (round == 11) {
				grapher.graph(t1, Vm1, Color.blue);
			}
			if (round == 12) {
				grapher.graph(t1, Vm1, Color.orange);
			}
			if (round == 13) {
				grapher.graph(t1, Vm1, Color.cyan);
			}
			if (round == 14) {
				grapher.graph(t1, Vm1, Color.magenta);
			}
			if (round == 15) {
				grapher.graph(t1, Vm1, Color.black);
			}
			if (round == 16) {
				grapher.graph(t1, Vm1, Color.yellow);
			}
			if (round == 17) {
				grapher.graph(t1, Vm1, Color.red);
			}
			if (round == 18) {
				grapher.graph(t1, Vm1, Color.gray);
			}
			if (round == 19) {
				grapher.graph(t1, Vm1, Color.lightGray);
			}
			if (round == 20) {
				grapher.graph(t1, Vm1, Color.green);
			}
			if (round == 21) {
				grapher.graph(t1, Vm1, Color.pink);
			}
			if (round == 22) {
				grapher.graph(t1, Vm1, Color.GREEN);
			}

			round++;
			grapher2.graph(time, arrayfreq, Color.GREEN);
		}
		grapher.display();
		grapher2.display();

	}
}