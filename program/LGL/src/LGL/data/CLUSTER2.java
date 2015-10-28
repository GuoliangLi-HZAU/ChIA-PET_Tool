/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

/**
 *
 * @author ligl
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Vector;

public class CLUSTER2 implements Comparable {

    private REGION2 head;
    private REGION2 tail;
    private double weight = 0.0;

    public CLUSTER2() {
        head = null;
        tail = null;
        weight = 0.0;
    }

    public CLUSTER2(PET pet, int extendLength) {
        if (pet.getHead().getStrand() == '+') {
            head = new REGION2(pet.getHead().getChrom(), pet.getHead().getLoci(), pet.getHead().getLoci() + extendLength, pet.getHead().getStrand());
        } else {
            head = new REGION2(pet.getHead().getChrom(), pet.getHead().getLoci() - extendLength, pet.getHead().getLoci(), pet.getHead().getStrand());
        }

        if (pet.getTail().getStrand() == '+') {
            tail = new REGION2(pet.getTail().getChrom(), pet.getTail().getLoci() - extendLength, pet.getTail().getLoci(), pet.getTail().getStrand());
        } else {
            tail = new REGION2(pet.getTail().getChrom(), pet.getTail().getLoci(), pet.getTail().getLoci() + extendLength, pet.getTail().getStrand());
        }
        sortHeadTail();
    }

    public CLUSTER2(REGION2 head, REGION2 tail) {
        this.head = head;
        this.tail = tail;
        this.weight = 1.0;
        sortHeadTail();
    }

    public CLUSTER2(REGION2 head, REGION2 tail, double weight) {
        this.head = head;
        this.tail = tail;
        this.weight = weight;
        sortHeadTail();
    }

    private void sortHeadTail() {
        if (head.compareTo(tail) > 0) {
            REGION2 tempRegion = head;
            head = tail;
            tail = tempRegion;
        }
    }

    @Override
    public String toString() {
        int sameChrom = 0;
        if (this.getHead().getChrom().compareTo(this.getTail().getChrom()) == 0) {
            sameChrom = 1;
        }
        return (this.getHead().toString() + "\t" + this.getTail().toString() + "\t" + sameChrom + "\t" + REGION2.calculateCenterDistance(this.getHead(), this.getTail()) + "\t" + this.getWeight());
    }

    public void combine(CLUSTER2 anotherCluster) throws Exception {
        this.head.combine(anotherCluster.getHead());
        this.tail.combine(anotherCluster.getTail());
        this.setWeight(this.getWeight() + anotherCluster.getWeight());
    }

    public int compareTo(Object anotherCluster) throws ClassCastException {
        if (!(anotherCluster instanceof CLUSTER2)) {
            throw new ClassCastException("A CLUSTER2 object expected.");
        }
        int result = this.getHead().compareTo(((CLUSTER2) anotherCluster).getHead());
        if (result == 0) {
            result = this.getTail().compareTo(((CLUSTER2) anotherCluster).getTail());
        }
        return result;
    }

    /**
     * @return the head
     */
    public REGION2 getHead() {
        return head;
    }

    /**
     * @param head the head to set
     */
    public void setHead(REGION2 head) {
        this.head = head;
    }

    /**
     * @return the tail
     */
    public REGION2 getTail() {
        return tail;
    }

    /**
     * @param tail the tail to set
     */
    public void setTail(REGION2 tail) {
        this.tail = tail;
    }

    public static Vector<CLUSTER2> loadCluster2(String Cluster2File) throws IOException {
        Vector<CLUSTER2> cluster2s = new Vector<CLUSTER2>();
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(Cluster2File))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String fields[] = line.split("\t");
            REGION2 head = new REGION2(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]), fields[3].charAt(0));
            REGION2 tail = new REGION2(fields[4], Integer.parseInt(fields[5]), Integer.parseInt(fields[6]), fields[7].charAt(0));

            double weight = 1.0;
            if (fields.length >= 9) {
                weight = Double.parseDouble(fields[8]);
            }

            CLUSTER2 cluster = new CLUSTER2(head, tail, weight);
            cluster2s.add(cluster);
        }
        fileIn.close();

        //Collections.sort(cluster2s);
        return cluster2s;
    }

    public static void saveCluster2(String Cluster2File, Vector<CLUSTER2> cluster2s) throws IOException {
        PrintWriter fileOut = new PrintWriter(new BufferedWriter(new FileWriter(Cluster2File)));

        for (int i = 0; i < cluster2s.size(); i++) {
            fileOut.println(cluster2s.elementAt(i).toString());
        }

        fileOut.close();
    }

    /**
     * @return the weight
     */
    public double getWeight() {
        return weight;
    }

    /**
     * @param weight the weight to set
     */
    public void setWeight(double weight) {
        this.weight = weight;
    }
}
