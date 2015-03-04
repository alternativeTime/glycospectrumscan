package unicarb;

import Amino_acids.glycan_value;
import Amino_acids.glycan_value_list;
import Amino_acids.glycans;
import GenerateImage.DatabaseMethods;
import org.apache.commons.fileupload.FileItem;
import org.apache.commons.io.IOUtils;
import org.junit.Test;

import javax.servlet.RequestDispatcher;
import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created by matthew on 04/03/15.
 */
public class GlycanMassCalculate {

    public static void printCombination(double[] nums, int[] counts,int count[], int
            startIndex, double totalAmont){


        if (startIndex == nums.length || totalAmont<=0.0001){
            //System.out.print("" + totalAmont + "=");
            int s=0;
            for(int j=0; j<count.length; j++){
                if(count[j]==0) s++;
            }
            if(s==8){
                for (int i=0; i<nums.length; i++){
                    //       System.out.print("" + counts[i] + "*" + nums[i] + " + ");

                    count[i]=counts[i];



                }
                //	 System.out.println("");
            }
            return ;
        }else  if (startIndex == nums.length -1){
            if (totalAmont%nums[startIndex]==0){
                counts[startIndex] = (int) (totalAmont/nums[startIndex]);
                //  printCombination(nums, counts,count, startIndex+1,0);

            }
        } else {
            for(int i =0; i<= (totalAmont/nums[startIndex]);i++){
                counts[startIndex] = i;
                printCombination(nums, counts,count, startIndex+1,totalAmont -
                        nums[startIndex]*i);
            }
        }


    }

    @Test
    public void loadFile() throws IOException {
        String output = "";
        InputStream reader = new FileInputStream(new File("/tmp/glycans.txt"));
        StringWriter writer = new StringWriter();
        IOUtils.copy(reader, writer);
        String theString = writer.toString();

        theString =theString.replaceAll("[\\t\\n\\r\\s]","-" );
        theString =theString.replaceAll("--","-" );

        ArrayList <Double> M_M = new ArrayList<Double>();
        ArrayList <Double> A_M = new ArrayList<Double>();

        String [] temp;
        String delimiter="-";
        temp = theString.split(delimiter);

        for(int i=0; i<temp.length; i++){
            if( (i+1) %2==0){
                A_M.add(     Double.parseDouble(temp[i]));
                //  System.out.println(i+": "+ Double.parseDouble(temp[i]));
            }else{
                M_M.add(Double.parseDouble(temp[i]));
            }
        }

        glycan_value_list gvl = new glycan_value_list();
        glycan_value gv_NeuGc = (glycan_value)gvl.glycan_value_list.get("NeuGc");
        glycan_value gv_NeuAc = (glycan_value)gvl.glycan_value_list.get("NeuAc");
        glycan_value gv_HexNAc = (glycan_value)gvl.glycan_value_list.get("HexNAc");
        glycan_value gv_Hexose = (glycan_value)gvl.glycan_value_list.get("Hexose");
        glycan_value gv_Deoxyhexose = (glycan_value)gvl.glycan_value_list.get("Deoxyhexose");
        glycan_value gv_Pentose = (glycan_value)gvl.glycan_value_list.get("Pentose");
        glycan_value gv_Sulphate = (glycan_value)gvl.glycan_value_list.get("Sulphate");
        glycan_value gv_Phosphate = (glycan_value)gvl.glycan_value_list.get("Phosphate");


        ArrayList <glycans> glycanss = new ArrayList();
        ArrayList <glycans> glycansNOInDB = new ArrayList();
        double nums[] = {79.9799,80.0642,132.1161,146.1430,162.1424,203.1950,291.2579,307.2573};
        String names[]={"Phosphate","Sulphate","Pentose","Deoxyhexose","Hexose","HexNAc","NeuAc","NeuGc"};
        int[] count = new int[nums.length];
        int[] counts = new int[nums.length];
        double totalnumbers[] = new double [A_M.size()];

        for(int i=0; i<totalnumbers.length; i++){
            totalnumbers[i] = A_M.get(i)+0.000001;
        }
        DatabaseMethods DBM = new DatabaseMethods();
        String [] CompositionInDB = DBM.GetComposition();


        for(int i=0; i<A_M.size(); i++){
            double result =0;

            int Phosphate = 0;
            int Sulphate = 0;
            int Pentose = 0;
            int Deoxyhexose = 0;
            int Hexose = 0;
            int HexNAc = 0;
            int NeuAc = 0;
            int NeuGc = 0;


            //
            // printCombination(nums, counts,count,0,totalnumbers[i]);
            String check = Integer.toString(count[0]) + Integer.toString(count[1])+Integer.toString(count[2])+Integer.toString(count[3])+Integer.toString(count[4])+Integer.toString(count[5])+Integer.toString(count[6])+Integer.toString(count[7]);

            int tc=0;
            /*for(int l=1; l< CompositionInDB.length; l++){

                if(CompositionInDB[l].equals(check) ){
                    tc ++ ;

                }
                if(tc==0 && l==CompositionInDB.length-1){

                    for(int n=0;n<count.length; n++){
                        System.out.print(count[n] +"-"+names[n]+" ");
                    }
                    System.out.println();

                    Phosphate = count[0];
                    Sulphate =count[1];
                    Pentose =count[2];
                    Deoxyhexose = count[3];
                    Hexose = count[4];
                    HexNAc = count[5];
                    NeuAc = count[6];
                    NeuGc = count[7];


                    glycans g = new glycans();
                    g.NeuGc = NeuGc;
                    g.NeuAc = NeuAc;
                    g.HexNAc = HexNAc;
                    g.Hexose = Hexose;
                    g.Deoxyhexose = Deoxyhexose;
                    g.Pentose = Pentose;
                    g.Sulphate = Sulphate;
                    g.Phosphate = Phosphate;
                    g.Monoisotopic_Mass = M_M.get(i);
                    g.Average_Mass = A_M.get(i);
                    glycansNOInDB.add(g);

                    String glycan_NoInDB =  Integer.toString (glycanss.size()+1) ;
                    request.setAttribute("glycan_NoInDB", glycan_NoInDB);
                    RequestDispatcher dispatcher = request.getRequestDispatcher("GlycanMassError.jsp");
                    if (dispatcher != null){

                        dispatcher.forward(request, response);


                    return;
                }
            }*/


            /*for(int l=1; l< CompositionInDB.length; l++){

                if(CompositionInDB[l].equals(check)){

                    Phosphate = count[0];
                    Sulphate =count[1];
                    Pentose =count[2];
                    Deoxyhexose = count[3];
                    Hexose = count[4];
                    HexNAc = count[5];
                    NeuAc = count[6];
                    NeuGc = count[7];


                    glycans g = new glycans();
                    g.NeuGc = NeuGc;
                    g.NeuAc = NeuAc;
                    g.HexNAc = HexNAc;
                    g.Hexose = Hexose;
                    g.Deoxyhexose = Deoxyhexose;
                    g.Pentose = Pentose;
                    g.Sulphate = Sulphate;
                    g.Phosphate = Phosphate;
                    g.Monoisotopic_Mass = M_M.get(i);
                    g.Average_Mass = A_M.get(i);
                    glycanss.add(g);

                }
            }*/


            for(int x=0; x<count.length; x++){
                count[x] =0;
            }

        }
        double total_monoisotopic_mass =0;
        double total_average_mass=0;
        for(glycans g: glycanss){
            total_monoisotopic_mass += g.Monoisotopic_Mass;
            total_average_mass +=g.Average_Mass;

        }
    } 


}
