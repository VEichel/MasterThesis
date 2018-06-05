package alternativeLIborModelTest;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.Iterator;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;

import alternativeLIborModelTest.ExcelData.DataSymbol;

public class AlternativeLiborModelWithExcelDataTest {

	public static void main(String[] args) {
		
		ExcelData[] excelData = readExcelData("C:\\Users\\vince\\Google Drive\\UNI\\MasterThesis\\Literatur\\Excel\\Market Data EUR USD.xls");
		
	}
	
	
	
	
	public static ExcelData[] readExcelData(String fileDirectory) {
		try {
	        FileInputStream file = new FileInputStream(new File(fileDirectory));
	
	        //Create Workbook instance holding reference to .xlsx file
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
}

