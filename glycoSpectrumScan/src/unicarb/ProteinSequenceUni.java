package unicarb;

import Amino_acids.NewMass_GS;
import Amino_acids.PeptideMass;
import Amino_acids.glycans;
import Amino_acids.glyco_sequence;
import glycoSpectrumScan.ProteinSequenceSplit;
import org.apache.commons.io.FileUtils;
import org.junit.Test;
import uk.ac.ebi.kraken.interfaces.uniprot.comments.Enzyme;
import uk.ac.ebi.kraken.uuw.services.remoting.EntryRetrievalService;
import uk.ac.ebi.kraken.uuw.services.remoting.UniProtJAPI;

import javax.servlet.RequestDispatcher;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;


/**
 * Servlet implementation class ProteinSequence_Servlet
 */

public class ProteinSequenceUni {
    private static final long serialVersionUID = 1L;
    glycoSpectrumScan.ProteinSequenceSplit PSS = new ProteinSequenceSplit();
    ArrayList<String> ProteinSequenceSplit;

    public ProteinSequenceUni() {
        super();

    }

    /*public int [] setMC(int Missed_cleavages){
        int [] MC;
        int [] MC_index = new int [Missed_cleavages+1];
        for(int i=0; i<MC_index.length; i++){
            MC_index[i] = ProteinSequenceSplit.size() - (Missed_cleavages-i);
        }
        ProteinSequenceSplit = PSS.missed_cleavages(ProteinSequenceSplit, Missed_cleavages);
        MC = new int [ProteinSequenceSplit.size()];
        MC = PSS.GetMC(MC_index, Missed_cleavages, MC);
        return MC;
    } */
    public static String getSequence(String accession){
        String sequence="";
        EntryRetrievalService entryRetrievalService = UniProtJAPI.factory.getEntryRetrievalService();
        //	UniProtEntry entry = (UniProtEntry) entryRetrievalService.getUniProtEntry(accession);
        //If entry with a given accession number is not found, entry will be equal null
        //	if (entry != null) {
        Object attribute = UniProtJAPI.factory.getEntryRetrievalService().getUniProtAttribute(accession ,  "ognl:sequence.value");
        sequence = attribute.toString();

        //	}

        return sequence;
    }

}

