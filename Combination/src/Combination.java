import java.util.LinkedList;


public class Combination {
	
	public static void main(String [ ] args)
	{

		double number =810.739; 
				
		double NeuGc=307.2573;
		double NeuAc=291.2579;
		double HexNAc=203.1950;
		double Hexose=162.1424;
		
		double Deoxyhexose=146.1430;
		double Pentose=132.1161;
		
		double Sulphate=80.0642;
		double Phosphate=79.9799;
		
		
		double count[] = {0,0,0,0,0,0,0,0};
		
		LinkedList <Double> Hexoses = new LinkedList <Double> ();
		LinkedList <Double> HexNAcs = new LinkedList <Double> ();
		LinkedList <Double> Deoxyhexoses = new LinkedList <Double> ();
		LinkedList <Double> Pentoses = new LinkedList <Double> ();
		LinkedList <Double> NeuAcs = new LinkedList <Double> ();
		LinkedList <Double> NeuGcs = new LinkedList <Double> ();
		LinkedList <Double> Phosphates = new LinkedList <Double> ();
		LinkedList <Double> Sulphates = new LinkedList <Double> ();
		
		double remainder=99999999.999999;
			for(;;){
				
				double number_of_NeuGc =  (int) (number / NeuGc);
				if(number_of_NeuGc!=0.0){
					Hexoses.add(number_of_NeuGc);	
				    remainder = number % NeuGc;
				}
				if(remainder == 0)break;
				
				//double number_of_NeuGc = (int)();
				
				
				
			//System.out.println(Hexoses_remainder);
			
			
			
			}
		
		
		
		
	}
}
