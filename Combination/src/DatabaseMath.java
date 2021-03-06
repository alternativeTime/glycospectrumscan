import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.LinkedList;


public class DatabaseMath {
	private static Connection getPostgreSQLConnection() throws Exception {
	 	Class.forName("org.postgresql.Driver");
	//    System.out.println("Driver Loaded.");
	 	String Schema_name="glycodb";
	 	String User_name="postgres";
	 	String Password="root";
	   
	 	return DriverManager.getConnection("jdbc:postgresql://localhost:5432/"+Schema_name,User_name,Password);
	}
	
	public static String[] GetComposition(){
		
		Connection conn=null;
		 try{
			 conn = getPostgreSQLConnection();
			 if(conn!=null){
				 String sql = "select phosphate, sulfate, pentose, deoxyhexose,hexose,hexnac,neuac,neugc  from composition where kdn=0 AND kdo=0 AND hexa=0 AND methyl=0 AND acetyl=0 AND other=0";
				 conn.setAutoCommit(false);  
				 PreparedStatement stmt = conn.prepareStatement(sql);
				 ResultSet rset = stmt.executeQuery();
				 LinkedList <Object[]> records=new LinkedList<Object[]>();
				 int cols = rset.getMetaData().getColumnCount();
				 while ( rset.next() )
			     {
							Object[] arr = new Object[cols];
						    for(int i=0; i<cols; i++){
						      arr[i] = rset.getObject(i+1);
						    }
						    records.add(arr);
						   
			      }
				 
				String [] results = new String [records.size()];
				for(int i=0; i<results.length; i++){
					results[i] ="";
				}
						 
				 for(int i=0; i<records.size();i++){
						for(int j=0; j<records.get(i).length; j++){
							results[i] += records.get(i)[j];
						//	CompositionID[i] = (String) records.get(i)[j];
						//	System.out.print(records.get(i)[j]+"  ");
						}
					//	System.out.println();
				} 
				 
				
				
				rset.close();
	    		stmt.close();
	    		conn.commit();
	    		conn.close();
	    		return results;
			 }
			 
		 }catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
		 }
		 return null;
	}

	public static void main(String[] args) {
		String [] results =GetComposition();
		for(String s: results){
			System.out.println(s);
		}
		
		
	}
}
