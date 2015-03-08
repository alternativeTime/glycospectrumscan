package test;

/**
 * Created by matthew on 07/03/15.
 */
import Amino_acids.*;
import org.apache.commons.io.FileUtils;
import org.junit.Test;
import uk.ac.ebi.kraken.uuw.services.remoting.EntryRetrievalService;
import uk.ac.ebi.kraken.uuw.services.remoting.UniProtJAPI;
import glycoSpectrumScan.ProteinSequenceSplit;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * Created by matthew on 6/03/15.
 */
public class glycanmassTest {


    @Test
    public void residue() throws IOException {

        ProteinSequenceSplit PSS = new ProteinSequenceSplit();
        String sequence = "";
        EntryRetrievalService entryRetrievalService = UniProtJAPI.factory.getEntryRetrievalService();

        Object attribute = UniProtJAPI.factory.getEntryRetrievalService().getUniProtAttribute("P01833", "ognl:sequence.value");
        sequence = attribute.toString();
        System.out.println(sequence);

        ArrayList <String> proteinSequenceSplit = new ArrayList();
        ArrayList <String> proteinSequenceSplit2 = new ArrayList();
        proteinSequenceSplit = PSS.ProteinSequence_Trypsin(sequence);

        PSS.missed_cleavages_calculate(proteinSequenceSplit, 2);


        int [] Hexoses = {3,5,3};
        int [] Deoxyhexoses = {1,0,0};
        int [] HexNAcs = {2,2,4};
        int [] Pentoses = {0,0,0};
        int [] NeuAcs  = {0,0,0};
        int [] Phosphates = {0,0,0};
        int [] NeuGcs  = {0,0,0};
        int [] Sulphates = {0,0,0};

        ArrayList<glycans> glycanss = new ArrayList();
        for(int i=0; i<Hexoses.length; i++){
            glycans g = new glycans();
            g.Hexose = Hexoses[i];
            g.HexNAc = HexNAcs[i];
            g.NeuAc = NeuAcs[i];
            g.NeuGc = NeuGcs[i];
            g.Deoxyhexose = Deoxyhexoses[i];
            g.Pentose = Pentoses[i];
            g.Phosphate = Phosphates[i];
            g.Sulphate = Sulphates[i];
            glycanss.add(g);
        }

        System.out.println("glycanss: " + glycanss.size());

        glycan_value_list gvl = new glycan_value_list();


        for(glycans g: glycanss){

            glycan_value gv_Hexose = (glycan_value)gvl.glycan_value_list.get("Hexose");
            double gv_Hexose_monoisotopic = gv_Hexose.monoisotopic;
            double gv_Hexose_average = gv_Hexose.average;

            glycan_value gv_HexNAc = (glycan_value)gvl.glycan_value_list.get("HexNAc");
            double gv_HexNAc_monoisotopic = gv_HexNAc.monoisotopic;
            double gv_HexNAc_average = gv_HexNAc.average;

            glycan_value gv_Deoxyhexose = (glycan_value)gvl.glycan_value_list.get("Deoxyhexose");
            double gv_Deoxyhexose_monoisotopic = gv_Deoxyhexose.monoisotopic;
            double gv_Deoxyhexose_average = gv_Deoxyhexose.average;

            glycan_value gv_Pentose = (glycan_value)gvl.glycan_value_list.get("Pentose");
            double gv_Pentose_monoisotopic = gv_Pentose.monoisotopic;
            double gv_Pentose_average = gv_Pentose.average;

            glycan_value gv_NeuAc = (glycan_value)gvl.glycan_value_list.get("NeuAc");
            double gv_NeuAc_monoisotopic = gv_NeuAc.monoisotopic;
            double gv_NeuAc_average = gv_NeuAc.average;

            glycan_value gv_NeuGc = (glycan_value)gvl.glycan_value_list.get("NeuGc");
            double gv_NeuGc_monoisotopic = gv_NeuGc.monoisotopic;
            double gv_NeuGc_average = gv_NeuGc.average;

            glycan_value gv_Phosphate = (glycan_value)gvl.glycan_value_list.get("Phosphate");
            double gv_Phosphate_monoisotopic = gv_Phosphate.monoisotopic;
            double gv_Phosphate_average = gv_Phosphate.average;

            glycan_value gv_Sulphate = (glycan_value)gvl.glycan_value_list.get("Sulphate");
            double gv_Sulphate_monoisotopic = gv_Sulphate.monoisotopic;
            double gv_Sulphate_average = gv_Sulphate.average;

            g.Average_Mass = g.Hexose * gv_Hexose_average + g.HexNAc * gv_HexNAc_average + g.Deoxyhexose * gv_Deoxyhexose_average + g.Pentose * gv_Pentose_average + g.NeuAc * gv_NeuAc_average + g.NeuGc * gv_NeuGc_average + g.Phosphate * gv_Phosphate_average + g.Sulphate * gv_Sulphate_average;
            g.Monoisotopic_Mass= g.Hexose * gv_Hexose_monoisotopic + g.HexNAc * gv_HexNAc_monoisotopic + g.Deoxyhexose * gv_Deoxyhexose_monoisotopic + g.Pentose * gv_Pentose_monoisotopic + g.NeuAc * gv_NeuAc_monoisotopic + g.NeuGc * gv_NeuGc_monoisotopic + g.Phosphate * gv_Phosphate_monoisotopic + g.Sulphate * gv_Sulphate_monoisotopic;
            System.out.println(g.Average_Mass);
            System.out.println(g.Monoisotopic_Mass);
        }

        double total_monoisotopic_mass =0;
        double total_average_mass=0;
        for(glycans g: glycanss){
            total_monoisotopic_mass += g.Monoisotopic_Mass;
            total_average_mass +=g.Average_Mass;
            System.out.println(g.Average_Mass);
            System.out.println(g.Monoisotopic_Mass);

        }

        //return glycanss;

        PeptideMass PM=PSS.getPeptideMass(proteinSequenceSplit, "M+H");

        glyco_sequence gsT = new glyco_sequence(PM.ProteinSequenceSplit,3,glycanss);

        ArrayList<String> findGs = gsT.FindGlycoSequence(PSS.missed_cleavages_calculate(proteinSequenceSplit, 1));
        //Combinations combinations = new Combinations();


        //
        //ArrayList<String> peptideMass = PSS.missed_cleavages_calculate(proteinSequenceSplit, 1);

        for(Iterator i = findGs.iterator();i.hasNext();) {
            //System.out.println(i.next());
            proteinSequenceSplit2.add(i.next().toString());
        }
        ProteinSequenceSplit PSS2 = new ProteinSequenceSplit();
        PeptideMass PM2=PSS2.getPeptideMass(proteinSequenceSplit2, "M+H");
        glyco_sequence gs = new glyco_sequence(PM2.ProteinSequenceSplit,3,glycanss);

        String Treat = "";
        Treat = "LODOACETAMID";
        if(Treat.equals("LODOACETAMID")){
			/* count how many time c appear in sequence*/

            for(int i=0; i<PM.ProteinSequenceSplit.size(); i++){
                int count = 0;
                int idx = 0;
                //		System.out.println(PM.ProteinSequenceSplit.get(i));
                //		System.out.print(PM.average_mass[i]+" ");
                //		System.out.println(PM.monoisotopic_mass[i]);

                while ((idx = PM.ProteinSequenceSplit.get(i).indexOf("C", idx)) != -1)
                {
                    idx++;
                    count++;
                }
                PM.average_mass[i] = PM.average_mass[i] + count*57.052;
                PM.monoisotopic_mass[i] = PM.monoisotopic_mass[i]+count*57.02146;

                //	System.out.println(PM.ProteinSequenceSplit.get(i));
                //	System.out.print(PM.average_mass[i]+" ");
                //	System.out.println(PM.monoisotopic_mass[i]);

            }
        }

        //from proteinsequence_serlet
        double [] glyco_sequence_m_mass = new double [gs.glyco_sequence.size()];
        double [] glyco_sequence_a_mass = new double [gs.glyco_sequence.size()];

        for(int i=0; i<gs.glyco_sequence.size();i++){
            for(int j=0; j<PM2.ProteinSequenceSplit.size();j++)
                if(gs.glyco_sequence.get(i).equals(PM2.ProteinSequenceSplit.get(j))){
                    glyco_sequence_m_mass[i] = PM2.monoisotopic_mass[j];
                    glyco_sequence_a_mass[i] = PM2.average_mass[j];
                    System.out.println("test mass print: " + PM2.ProteinSequenceSplit.get(j) + " " + PM2.monoisotopic_mass[j] );
                }
        }

        System.out.println("check: " + glyco_sequence_m_mass[1]);
        LinkedList <NewMass_GS> gs_NM = new LinkedList();

        for(int i=0; i<gs.Combination.size();i++){
            //	System.out.println("Mass before adding to " +gs.glyco_sequence.get(i)+" is:"+ glyco_sequence_m_mass[i]+" "+ glyco_sequence_a_mass[i]);

            NewMass_GS ngs = new NewMass_GS();

            for(int j=0; j<gs.Combination.get(i).Combination.size(); j++){
                ngs.new_m_mass.add(glyco_sequence_m_mass[i]+gs.Combination.get(i).m_mass[j]);
                ngs.new_a_mass.add(glyco_sequence_a_mass[i]+gs.Combination.get(i).a_mass[j]);

                //		System.out.println(gs.Combination.get(i).Combination.get(j)+" : " +gs.Combination.get(i).m_mass[j]+" "+ gs.Combination.get(i).a_mass[j]);
            }
            gs_NM.add(ngs);
        }



        String write_file="";

        for(int i=0; i<gs_NM.size(); i++){
            write_file += "Mass before adding to " +gs.glyco_sequence.get(i)+" is:"+ glyco_sequence_m_mass[i]+" "+ glyco_sequence_a_mass[i]+"\n";
            //	System.out.println("Mass before adding to " +gs.glyco_sequence.get(i)+" is:"+ glyco_sequence_m_mass[i]+" "+ glyco_sequence_a_mass[i]);

            for(int j=0; j<gs_NM.get(i).new_m_mass.size(); j++){
                write_file += gs.Combination.get(i).Combination.get(j)+": New Monoisotopic mass:"+gs_NM.get(i).new_m_mass.get(j)+" New Average mass:"+gs_NM.get(i).new_a_mass.get(j)+"\n";
                //	System.out.println(gs.Combination.get(i).Combination.get(j)+" : " +gs.Combination.get(i).m_mass[j]+" "+ gs.Combination.get(i).a_mass[j]);
                //	System.out.println("new m mass is : "+ gs_NM.get(i).new_m_mass.get(j)+" new a mass is: "+ gs_NM.get(i).new_a_mass.get(j)  );
            }

        }


        String tempFile = "/tmp/Result.txt";
        //Delete if tempFile exists
        File fileTemp = new File(tempFile);
        if (fileTemp.exists()){
            fileTemp.delete();
        }

        FileUtils.writeStringToFile(new File("/tmp/Result.txt"), write_file);





    }

}

