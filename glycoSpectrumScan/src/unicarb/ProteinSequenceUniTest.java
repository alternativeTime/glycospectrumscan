package unicarb;

import Amino_acids.*;
import glycoSpectrumScan.ProteinSequenceSplit;
import org.apache.commons.logging.Log;
import org.junit.Test;
import uk.ac.ebi.kraken.uuw.services.remoting.EntryRetrievalService;
import uk.ac.ebi.kraken.uuw.services.remoting.UniProtJAPI;

import java.util.ArrayList;
import java.util.logging.Logger;

import static org.junit.Assert.*;

public class ProteinSequenceUniTest {

    ProteinSequenceSplit PSS = new ProteinSequenceSplit();
    ArrayList<String> ProteinSequenceSplit;
    private int Missed_cleavages = 2;

    @Test
    public void getSequence(){
        String sequence="";
        String accession = "P01833";
        EntryRetrievalService entryRetrievalService = UniProtJAPI.factory.getEntryRetrievalService();
        //	UniProtEntry entry = (UniProtEntry) entryRetrievalService.getUniProtEntry(accession);
        //If entry with a given accession number is not found, entry will be equal null
        //	if (entry != null) {
        Object attribute = UniProtJAPI.factory.getEntryRetrievalService().getUniProtAttribute(accession ,  "ognl:sequence.value");
        sequence = attribute.toString();
        System.out.println("sequence: " + sequence);

        String Enzyme = "TRYPSIN";
        if(Enzyme.equals("TRYPSIN")){
            ProteinSequenceSplit = PSS.ProteinSequence_Trypsin(sequence);
            if(ProteinSequenceSplit.size() <= Missed_cleavages){
                Missed_cleavages = ProteinSequenceSplit.size()-1;
            }

            //MC = setMC(Missed_cleavages);
            //request.setAttribute("MC",MC);
        }

        for(String p : ProteinSequenceSplit) {
            System.out.println("\n" + p.toString());

        }

        PeptideMass PM= PSS.getPeptideMass(ProteinSequenceSplit, "M+H");
        String Treat = "LODOACETAMID";
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


        //need to create glycan masses see glycanmass_serlets
        //glyco_sequence to calculate options using peptidelist etc

        //return sequence;
    }

}