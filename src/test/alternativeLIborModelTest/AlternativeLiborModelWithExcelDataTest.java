package alternativeLIborModelTest;

import alternativeLIborModelTest.ExcelData.DataSymbol;
import net.finmath.stochastic.RandomVariableInterface;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;

import java.io.*;

public class AlternativeLiborModelWithExcelDataTest {

	public static void main(String[] args) {
		
		ExcelData[] excelData = readExcelData("C:\\Users\\vince\\Google Drive\\UNI\\MasterThesis\\Literatur\\Excel\\Market Data EUR USD.xls");
	}
	
	
	
	
	public static ExcelData[] readExcelData(String fileDirectory) {
		try {
	        FileInputStream file = new FileInputStream(new File(fileDirectory));
	
	        //Create Workbook instance holding reference to .xls file
	        HSSFWorkbook workbook = new HSSFWorkbook(file);

	        //Get first/desired sheet from the workbook
	        HSSFSheet sheet = workbook.getSheetAt(0);

	        //Iterate through each rows one by one
	        int rows = sheet.getPhysicalNumberOfRows();
	        ExcelData[] excelData = new ExcelData[rows-3];
	       
	        for (int rowIndex = 3; rowIndex < rows; rowIndex++) {
	            Row row = sheet.getRow(rowIndex);
	           
	            
	            excelData[rowIndex - 3] = new ExcelData(row.getCell(0).getDateCellValue(), 
	            									DataSymbol.valueOf(row.getCell(1).getStringCellValue()),
	            									row.getCell(2).getNumericCellValue());
	        }
	        workbook.close();
	        file.close();
	        return excelData;
	       
		} catch (FileNotFoundException fileNotFoundException) {
			System.out.println("The Excel Data File was not found! Change the readExcelData input to match your filepath.");
			fileNotFoundException.printStackTrace();
	    } catch (Exception e) {
	        e.printStackTrace();
	    }
		return null;
	}
	
	public static void printDateToExcel(RandomVariableInterface rv, int colum, String fileName) {
		
		String fileDirectory = "C:\\Users\\vince\\Desktop" + "\\" + fileName + ".xls";
		
		
		try {
			File file    = new File(fileDirectory);
			HSSFWorkbook workbook;
			HSSFSheet sheet;
			String sheetname = "test";
			boolean oldFile = false; 
			
			if(file.exists()) {
				FileInputStream fileInput = new FileInputStream(file);
				
				workbook = new HSSFWorkbook(fileInput);
				sheet = workbook.getSheetAt(0);
				fileInput.close();
				oldFile = true;
			} else {
				workbook = new HSSFWorkbook();
				sheet = workbook.createSheet(sheetname);
			}	         
   
	        for (int rowNumber = 0; rowNumber < rv.size(); rowNumber++) {
	        	Row row = oldFile ? sheet.getRow(rowNumber) : sheet.createRow(rowNumber);
	        	Cell cell = row.createCell(colum);
	        	cell.setCellValue(rv.get(rowNumber));
	        }
	        FileOutputStream fileOut = new FileOutputStream(file);
	        workbook.write(fileOut);
	        workbook.close();
	        fileOut.close();
	        
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}

