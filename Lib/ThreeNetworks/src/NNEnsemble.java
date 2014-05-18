import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;

import org.encog.neural.activation.ActivationSigmoid;
import org.encog.neural.activation.ActivationTANH;
import org.encog.neural.data.NeuralData;
import org.encog.neural.data.NeuralDataPair;
import org.encog.neural.data.NeuralDataSet;
import org.encog.neural.data.basic.BasicNeuralDataSet;
import org.encog.neural.networks.BasicNetwork;
import org.encog.neural.networks.layers.BasicLayer;
import org.encog.neural.networks.logic.FeedforwardLogic;
import org.encog.neural.networks.training.Train;
import org.encog.neural.networks.training.lma.LevenbergMarquardtTraining;
import org.encog.neural.networks.training.propagation.back.Backpropagation;
import org.encog.neural.networks.training.strategy.RequiredImprovementStrategy;
import org.encog.util.logging.Logging;

public class NNEnsemble {
//	private double trainInput[][];
//	private double trainIdeal[][];
//	
//	private double testInput[][];
//	private double testIdeal[][];
	
	private String trainFile, testFile, outputFile;
	private NeuralDataSet trainingSet, testingSet;
	
	private double MAX_ERROR = 1e-3; 
	private int inNum = 4;
	public NNEnsemble(int inNum, String train, String test, String output) {
		this.inNum = inNum;
		this.trainFile = train;
		this.testFile = test;
		this.outputFile = output;
	}
	public NNEnsemble() {
		this(4, "test/data0.txt", "test/test0.txt", "test/output0.txt");
	}
	private NeuralDataSet readInData(String fileName) {
		try {
			// Count number of lines
			LineNumberReader reader = new LineNumberReader(new FileReader(fileName));
			reader.skip(Long.MAX_VALUE);
			reader.close();
			int lineNumber = reader.getLineNumber();
			
			// Fill in training data
			reader = new LineNumberReader(new FileReader(fileName));
			String line;
			String[] lines;
			System.out.println(lineNumber);
			double[][] input = new double[lineNumber][inNum];
			double[][] ideal = new double[lineNumber][1];
			
			for (int i = 0; i < lineNumber; i++) {
				line = reader.readLine();
				lines = line.split("[: ]");
				// 0(target value) 2 4 6 8
				for (int j = 0; j < inNum; j++)
					input[i][j] = Double.parseDouble(lines[(j + 1) * 2]);
				ideal[i][0] = Double.parseDouble(lines[0]);
			}
			reader.close();
			
			NeuralDataSet dataSet = new BasicNeuralDataSet(input, ideal);
			return dataSet;
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
	}
	private void printDouble(double[][] data) {
		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data[0].length; j++)
				System.out.print(data[i][j] + " ");
			System.out.println();
		}
	}
	public void run() throws IOException {
		trainingSet = readInData(trainFile);
		testingSet = readInData(testFile);
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
		Logging.stopConsoleLogging();
		
		BasicNetwork network = new BasicNetwork();
		// does threshold affects?
		network.addLayer(new BasicLayer(new ActivationSigmoid(),true,inNum));
		network.addLayer(new BasicLayer(new ActivationSigmoid(),true,inNum));
		network.addLayer(new BasicLayer(new ActivationSigmoid(),true,1));
		network.setLogic(new FeedforwardLogic());
		network.getStructure().finalizeStructure();
		network.reset();

		// train the neural network
		final Train train = new LevenbergMarquardtTraining(network, trainingSet);
		
		// reset if improve is less than 1% over 5 cycles
		train.addStrategy(new RequiredImprovementStrategy(5));

		int epoch = 1;

		double prev;
		train.setError(Double.MAX_VALUE);
		
		// simple stopping case with one minimal (without momentum) 
		do {
			prev = train.getError();
			train.iteration();
			System.out
					.println("Epoch #" + epoch +
							" Error:" + train.getError());
			epoch++;
		} while(train.getError() > MAX_ERROR
				&& prev > train.getError());

		// compute the results for the test set
//		System.out.println("Neural Network Results:");
//		for(NeuralDataPair pair: testingSet ) {
//			final NeuralData output = network.compute(pair.getInput());
//			System.out.println("actual=" + output.getData(0) + ",ideal=" + pair.getIdeal().getData(0));
//		}
		for(NeuralDataPair pair: testingSet ) {
			final NeuralData output = network.compute(pair.getInput());
			bw.write(output.getData(0) + "\n");
		}
		bw.close();

	}
	// Usage:
	// NNEnsemble [train] [test] [output]
	public static void main(String[] args) {
		NNEnsemble neuralEn;
		neuralEn = args.length < 4 ? new NNEnsemble(4, args[0], args[1], args[2]) 
		: new NNEnsemble();
		try {
			neuralEn.run();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
