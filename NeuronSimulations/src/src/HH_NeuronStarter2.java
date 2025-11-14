package src;


// Updated by Rebecca Mantione, orginal code from Chris Fietkiewicz. Hodgkin-Huxley model of a neuron.
// NOTE: Requires SwingGraphics.java for graphing.
import java.awt.*;
import java.util.ArrayList;
import java.util.Scanner;

public class HH_NeuronStarter2 {
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
	static int found = 0;
	static double save = 0;
	static int first = 0;

// Constructor
	public HH_NeuronStarter2(double dt) {
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

	public static void simulate(double dt, int N, double[] t, double[] Vm, double stimAmplitude, double stimDuration) {
		HH_NeuronStarter n = new HH_NeuronStarter(dt); // Create a neuron
		double firstMsec = 0.5; // Time for 1st stimulus
		int first = (int) Math.round(firstMsec / n.getdt());
		int intStimDuration = (int) Math.round(stimDuration / n.getdt()) - 1;
//		System.out.println("Simulation settings: amplitude = " + stimAmplitude + ", dt = " + dt + " msec...");
		for (int i = 0; i < N; i++) {
			t[i] = i * n.getdt();
			if (i >= first && i <= (first + intStimDuration)) {
				Vm[i] = n.calculateNextTimeStep(stimAmplitude);
//				System.out.println("Found an action potential!");
			}else // Normally, recalcuate membrane potential with no stimulus
				Vm[i] = n.calculateNextTimeStep(0.0);
		}
	}

// Main method that runs the simulation
	public static void main(String[] args) {
		ArrayList<double[]> arrayList = new ArrayList<>();
		Scanner reader = new Scanner(System.in); // Reading from System.in
		System.out.println("Enter starting amplitude: ");
		double n = reader.nextInt(); // Scans the next token of the input as an int
		// Once finished
		reader.close();
//		System.out.println("readin" + n);
		double dt = 0.01; // Time step duration
		double duration = 10.0; // Simulation duration
		double stimDuration = 0.1; // Stimulus duration
// Experiment
		double stimAmplitude = n;
		int N = (int) (duration / dt);
		int round = 0;
		SwingGraphics grapher = new SwingGraphics(); // Begin graphing
		
		for (double i = n; i < n+11; i++) {
			stimAmplitude = i;
			System.out.println("Simulation Settings: amplitude = " + stimAmplitude);
			double[] t1 = new double[N];
			double[] Vm1 = new double[N];
			simulate(dt, N, t1, Vm1, stimAmplitude, stimDuration);
			for(int k=0; k< Vm1.length; k++) {
			if (Vm1[k] >= 0) {
				found++;
				first++;
				if (found == 1) {
					System.out.println("Found action potential!");
				}
				if (first == 1) {
					save = stimAmplitude;
				}
			}
			}
			found = 0;
			if(round == 1) {
				grapher.graph(t1, Vm1, Color.GREEN);
			}
			else if(round == 2) {
				grapher.graph(t1, Vm1, Color.PINK);
			}
			else if(round == 3) {
				grapher.graph(t1, Vm1, Color.ORANGE);
			}
			else if(round == 4) {
				grapher.graph(t1, Vm1, Color.CYAN);
			}
			else if(round == 5) {
				grapher.graph(t1, Vm1, Color.MAGENTA);
			}
			else if(round == 6) {
				grapher.graph(t1, Vm1, Color.BLACK);
			}
			else if(round == 7) {
				grapher.graph(t1, Vm1, Color.YELLOW);
			}
			else if(round == 8) {
				grapher.graph(t1, Vm1, Color.RED);
			}
			else if(round == 9) {
				grapher.graph(t1, Vm1, Color.GRAY);
			}
			else if (round == 10) {
				grapher.graph(t1, Vm1, Color.LIGHT_GRAY);
			}
			round ++;

		}
		
		
		double stimAmplitude2 = 66;
		int N2 = (int) (duration / dt);
		double[] t1 = new double[N];
		double[] Vm1 = new double[N];
		simulate(dt, N2, t1, Vm1, stimAmplitude2, stimDuration);
//		grapher.graph(t1, Vm1, Color.GREEN);
		grapher.display();
		if (save == 0.0 ) {
			System.out.println("No action potential found");
		}else {
		System.out.println("Input threshold was " + save);
		}
	}
}
