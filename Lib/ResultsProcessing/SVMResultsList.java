package ResultsProcessing;

import java.io.*;

public class SVMResultsList {
	public static void main(String args[]) {
		BufferedReader br, br2;
		BufferedWriter bw;
		try {
			br = new BufferedReader(new FileReader("record.txt"));
			bw = new BufferedWriter(new FileWriter("result.txt", false));
			for (int i = 0; i < 10; i++) {
				br2 = new BufferedReader(new FileReader("result" + i + ".txt"));
				String s, s2;
				while ((s = br2.readLine()) != null) {
					s2 = br.readLine();
					bw.write(s + " " + s2);
					bw.newLine();
				}
				br2.close();
			}
			
			br.close();
			bw.close();
		} catch (IOException e) {
			System.out.println("record.txt or result[i].txt Not Found. Or other Reasons.");
		}
	}
}
