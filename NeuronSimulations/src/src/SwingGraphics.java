package src;

// Written by Ian Waldschmidt at Case Western Reserve University
// and Chris Fietkiewicz at Hobart and William Smith Colleges. The "graph" method
// receives arrays for x-axis and y-axis values, plus a line color. You can graph
// multilple lines and then call the "display" method once to display all of the lines.
// A main() method is provided for demonstration, but the methods are public and can
// be used in other programs by declaring an instance of the Grapher class and calling
// the graph() and display() methods in the other program.
// Notes: (1) "graph" must be called for all data before "display".
//        (2) import java.awt.* to use the Color class in other programs.

import java.util.*;
import java.awt.*;
import java.awt.font.*;
import java.awt.geom.*;
import javax.swing.*;
import java.text.*; // For DecimalFormat

public class SwingGraphics extends JPanel {
    JFrame f = new JFrame();
    private ArrayList<double[]> xdatas = new ArrayList<double[]>();
    private ArrayList<double[]> ydatas = new ArrayList<double[]>();
    private ArrayList<Color> colors = new ArrayList<Color>();

    final int PAD = 30;

    public void graph(double[] xdata, double[] ydata, Color color){
        if(xdata.length != ydata.length){
            System.out.println("Input arrays must be the same length");
            return;
        }

        this.xdatas.add(xdata);
        this.ydatas.add(ydata);
        this.colors.add(color);
    }

    public void display(){
        f.add(this);
        f.setSize(1000,600);
        f.setLocation(100,100);
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); // Ensure that program terminates when window is closed
        f.setVisible(true);
    }

    protected void paintComponent(Graphics g){
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D)g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        int w = getWidth();
        int h = getHeight();
        //Draw ordinate
        g2.draw(new Line2D.Double(PAD, PAD, PAD, h-PAD));
        //Draw abcissa
        g2.draw(new Line2D.Double(PAD, h-PAD, w-PAD, h-PAD));
        //Draw labels
        Font font = g2.getFont();
        FontRenderContext frc = g2.getFontRenderContext();
        LineMetrics lm = font.getLineMetrics("0", frc);
        float sh = lm.getAscent() + lm.getDescent();

        //Draw data
        double xMin = getMin(xdatas);
        double yMin = getMin(ydatas);
        double xRange = getMax(xdatas) - xMin;
        double yRange = getMax(ydatas) - yMin;
        double xScale = (double)(w - 2*PAD) / Math.abs(xRange);
        double yScale = (double)(h - 2*PAD) / Math.abs(yRange);
   		//System.out.println(xMin + ", " + yMin + ", " + xRange + ", " + yRange + ", " + xScale + ", " + yScale);
        for(int j = 0; j < xdatas.size(); j++){
            double[] xdata = this.xdatas.get(j);
            double[] ydata = this.ydatas.get(j);
            Color color = this.colors.get(j);
            //Draw lines for each data set
            g2.setPaint(color);
            for(int i = 0; i < xdata.length-1; i++) {
                double x1 = PAD + xScale * (xdata[i] - xMin);
                double y1 = h - PAD - yScale * (ydata[i] - yMin);
                double x2 = PAD + xScale * (xdata[i+1] - xMin);
                double y2 = h - PAD - yScale * (ydata[i+1] - yMin);
                g2.draw(new Line2D.Double(x1, y1, x2, y2));
            }
            //Draw points for each data set
            g2.setPaint(color);
            for(int i = 0; i < xdata.length; i++) {
                double x = PAD + xScale * (xdata[i] - xMin);
                double y = h - PAD - yScale * (ydata[i] - yMin);
                g2.fill(new Ellipse2D.Double(x-2, y-2, 4, 4));
            }
        }

        //Write x axis values
        g2.setPaint(Color.BLACK);
        double[] xaxis = new double[5];
        xaxis[0] = getMin(xdatas);
        xaxis[4] = getMax(xdatas);
        xaxis[1] = xRange / 4 + xaxis[0];
        xaxis[2] = xRange / 2 + xaxis[0];
        xaxis[3] = xRange * 3 / 4 + xaxis[0];
        int xcoord = PAD - 4;
        int ycoord = h - (PAD/4);
        for(int i = 0; i < 5; i++){
			DecimalFormat form = new DecimalFormat("0.#");
            String value = form.format(xaxis[i]) + "";
            g2.drawString(value, xcoord, ycoord);
            xcoord += (w - PAD) / 4 - PAD + 22;
        }

        //Write y axis values
        double[] yaxis = new double[5];
        yaxis[0] = getMin(ydatas);
        yaxis[4] = getMax(ydatas);
        yaxis[1] = yRange / 4 + yaxis[0];
        yaxis[2] = yRange / 2 + yaxis[0];
        yaxis[3] = yRange * 3 / 4 + yaxis[0];
        xcoord = PAD/4;
        ycoord = h - PAD + 4;
        for(int i = 0; i < 5; i++){
			DecimalFormat form = new DecimalFormat("0.#");
            String value = form.format(yaxis[i]) + "";
            g2.drawString(value, xcoord, ycoord);
            ycoord -= (h - PAD) / 4 - PAD + 22;
        }
    }

    private double getMin(ArrayList<double[]> data){
        double min = Integer.MAX_VALUE;
        for(int j = 0; j < data.size(); j++){
            for(int i = 0; i < data.get(j).length; i++){
                if(data.get(j)[i] < min){
                    min = data.get(j)[i];
                }
            }
        }
        return min;
    }

    private double getMax(ArrayList<double[]> data){
        double max = -Integer.MAX_VALUE;
        for(int j = 0; j < data.size(); j++){
            for(int i = 0; i < data.get(j).length; i++){
                if(data.get(j)[i] > max){
                    max = data.get(j)[i];
                }
            }
        }
        return max;
    }

    public static void main(String[] args){
        SwingGraphics grapher = new SwingGraphics(); // Create single Grapher for all lines

        // Make square with vertices at (2,2), (2,3), (3,3) and (3, 2)
        double[] xSquare = {2, 2, 3, 3, 2}; // x-coordinates for all sides of a square
        double[] ySquare = {2, 3, 3, 2, 2}; // y-coordinates for all sides of a square

        // Make triangle with vertices at (0,1), (1,2), and (1, 0)
        double[] xTriangleSide1 = {0, 1}; // x-coordinates for side #1 of a triangle
        double[] yTriangleSide1 = {1, 2}; // y-coordinates for side #1 of a triangle
        double[] xTriangleSide2 = {1, 1}; // x-coordinates for side #2 of a triangle
        double[] yTriangleSide2 = {2, 0}; // y-coordinates for side #2 of a triangle
        double[] xTriangleSide3 = {1, 0}; // x-coordinates for side #3 of a triangle
        double[] yTriangleSide3 = {0, 1}; // y-coordinates for side #3 of a triangle


        // Prepare all graphics
        grapher.graph(xSquare, ySquare, Color.BLUE);
        grapher.graph(xTriangleSide1, yTriangleSide1, Color.BLUE);
        grapher.graph(xTriangleSide2, yTriangleSide2, Color.RED);
        grapher.graph(xTriangleSide3, yTriangleSide3, Color.GREEN);

        grapher.display(); // Display all lines when ready
    }
}
