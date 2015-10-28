/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Calendar;
import java.util.Collections;
import java.util.Vector;

/**
 *
 * @author ligl
 */
public class PET implements Comparable {

    private HIT2 head;
    private HIT2 tail;
    private double weight = 0.0;
    private int index;

    public PET() {
        head = null;
        tail = null;
        weight = 0.0;
        index = -1;
    }

    public PET(HIT2 hit1, HIT2 hit2, double weight) {
        this.head = hit1;
        this.tail = hit2;
        this.weight = weight;
        this.index = -1;
    }

    //using in MergeSimilarPETs
    public PET(HIT2 hit1, HIT2 hit2, double weight, int index, int sortLabel) {
        if (sortLabel == 0) {
            this.head = hit1;
            this.tail = hit2;
        } else if (sortLabel > 0) {
            if (hit1.compareTo(hit2) <= 0) {
                this.head = hit1;
                this.tail = hit2;
            } else {
                this.head = hit2;
                this.tail = hit1;
            }
        } else { // sortLabel < 0
            if (hit1.compareTo(hit2) > 0) {
                this.head = hit1;
                this.tail = hit2;
            } else {
                this.head = hit2;
                this.tail = hit1;
            }
        }
        this.weight = weight;
        this.index = -1;
    }

    public PET(HIT2 hit1, HIT2 hit2, double weight, int index) {
        this.head = hit1;
        this.tail = hit2;
        this.weight = weight;
        this.index = index;
    }

    @Override
    public String toString() {
        return (this.getHead().toString() + "\t" + this.getTail().toString() + "\t" + this.getWeight() + "\t" + this.getIndex());
    }

    public int compareTo(Object anotherPet) throws ClassCastException {
        if (!(anotherPet instanceof PET)) {
            throw new ClassCastException("A PET object expected.");
        }
        int result = getHead().compareTo(((PET) anotherPet).getHead());
        if (result == 0) {
            result = getTail().compareTo(((PET) anotherPet).getTail());
        }
        return result;
    }

    /**
     * @return the head
     */
    public HIT2 getHead() {
        return head;
    }

    /**
     * @param head the head to set
     */
    public void setHead(HIT2 head) {
        this.head = head;
    }

    /**
     * @return the tail
     */
    public HIT2 getTail() {
        return tail;
    }

    /**
     * @param tail the tail to set
     */
    public void setTail(HIT2 tail) {
        this.tail = tail;
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

    public int getSpan() {
        int distance = Integer.MAX_VALUE;
        if (this.getHead().getChrom().compareTo(this.getTail().getChrom()) == 0) // the same chromosome
        {
            if (this.getHead().getStrand() == this.getTail().getStrand()) // same strand
            {
                distance = Math.abs(this.getHead().getLoci() - this.getTail().getLoci());
            }
        }
        return distance;
    }

    public static Vector<PET> load(String petFile) throws IOException {
        System.out.println("start loading PETs from " + petFile);
        Vector<PET> pets = new Vector<PET>();
        Calendar rightNow = Calendar.getInstance();

        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(petFile))));
        String line;
        int nPETs = 0;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String[] fields = line.split("\t");
            HIT2 hit1 = new HIT2(fields[0], Integer.parseInt(fields[1]), fields[2].charAt(0));
            HIT2 hit2 = new HIT2(fields[3], Integer.parseInt(fields[4]), fields[5].charAt(0));
            double weight = 1.0;
            if (fields.length >= 7) {
                weight = Double.parseDouble(fields[6]);
            }
            PET newPet = new PET(hit1, hit2, weight);
            if (fields.length >= 8) {
                newPet.setIndex(Integer.parseInt(fields[7]));
            }
            pets.add(newPet);
            nPETs++;
            if (nPETs % 1000000 == 0) {
                rightNow = Calendar.getInstance();
                System.out.println("[" + rightNow.getTime().toString() + "] " + (nPETs / 1000000) + "M PETs read from " + petFile);
            }
        }
        fileIn.close();
        System.out.println("Totally, " + nPETs + " PETs read from " + petFile);

        return pets;
    }

    public static void save(String petFile, Vector<PET> pets) throws IOException {
        PrintWriter fileOut = new PrintWriter(new BufferedWriter(new FileWriter(petFile)));

        for (int i = 0; i < pets.size(); i++) {
            fileOut.println(pets.elementAt(i).toString());
            if (i % 1000000 == 0) {
                fileOut.flush();
            }
        }
        fileOut.close();
    }

    public static void save(String petFile, Vector<PET> pets, boolean append) throws IOException {
        PrintWriter fileOut = new PrintWriter(new BufferedWriter(new FileWriter(petFile, append)));

        for (int i = 0; i < pets.size(); i++) {
            fileOut.println(pets.elementAt(i).toString());
            if (i % 1000000 == 0) {
                fileOut.flush();
            }
        }
        fileOut.close();
    }

    public static void addPETs(Vector<PET> pets, PairedTagMaps pairedTagMaps) {
        Vector<MAP> bestHeadMaps = pairedTagMaps.getHeadMaps().bestMaps();
        Vector<MAP> bestTailMaps = pairedTagMaps.getTailMaps().bestMaps();
        if ((bestHeadMaps.size() <= 0) || (bestTailMaps.size() <= 0)) { // return if either head or tail does not have mappings
            return;
        }
        //double weight = 1.0 * pairedTagMaps.getCount() / (bestHeadMaps.size() * bestTailMaps.size());
        double weight = 1.0 / (bestHeadMaps.size() * bestTailMaps.size());

        for (int i = 0; i < bestHeadMaps.size(); i++) {
            for (int j = 0; j < bestTailMaps.size(); j++) {
                addPET(pets, bestHeadMaps.elementAt(i).toHIT2(), bestTailMaps.elementAt(j).toHIT2(), weight);
            }
        }
    }

    public static void addPET(Vector<PET> pets, HIT2 hit1, HIT2 hit2, double weight) {
        if ((hit1 == null) && (hit2 == null)) {
            return;
        }
        PET pet = new PET(hit1, hit2, weight);
        pets.add(pet);
    }

    // The purpose of the sort function here is for removing the PCR redundancy
    // mainly used in MapWeight class
    public static void addPET(Vector<PET> pets, HIT2 hit1, HIT2 hit2, double weight, int index, int sortLabel) {
        if ((hit1 == null) && (hit2 == null)) {
            return;
        }
        PET pet = new PET(hit1, hit2, weight, index, sortLabel);
        pets.add(pet);
    }

    public static Vector<PET> mergeSamePets(Vector<PET> pets) {
        if (pets.size() <= 1) {
            return pets;
        }

        Vector<PET> newChiapets = new Vector<PET>(); // Interaction-PET
        Collections.sort(pets);
        PET pet = pets.elementAt(0);
        for (int i = 1; i < pets.size(); i++) {
            if (pet.compareTo(pets.elementAt(i)) == 0) {
                pet.setWeight(pet.getWeight() + pets.elementAt(i).getWeight());
            } else {
                newChiapets.add(pet);
                pet = pets.elementAt(i);
            }
        }
        newChiapets.add(pet);
        return newChiapets;
    }

    // difference from previous version
    // 1. the previous version just works with the continuous coordinate
    //    can't process the following example
    //    chr1  1   +  chr1  100 +
    //    chr1  1   +  chr1  1000000 +
    //    chr1  2   +  chr1  100 +
    public static Vector<PET> mergeSimilarPets(Vector<PET> pets, int distance) {
        if (pets.size() <= 1) {
            return pets;
        }

        Vector<PET> newChiapets = new Vector<PET>(); // Interaction-PET
        Collections.sort(pets);
        while (pets.size() >= 1) {
            PET pet = pets.elementAt(0);
            pets.remove(0);
            for (int i = 0; i < pets.size();) {
                if (HIT2.calculateDistance(pet.getHead(), pets.elementAt(i).getHead()) <= distance) {
                    if (HIT2.calculateDistance(pet.getTail(), pets.elementAt(i).getTail()) <= distance) {
                        double weight = pet.getWeight() + pets.elementAt(i).getWeight();
                        pet = pets.elementAt(i);
                        pet.setWeight(weight);
                        pets.remove(i);
                    } else {
                        i++;
                    }
                } else {
                    break;
                }
            }
            newChiapets.add(pet);
        }
        return newChiapets;
    }//*/

    
    public static Vector<Integer> PetsDistnaces(Vector<PET> pets) {
        Vector<Integer> distances = new Vector<Integer>();
        if (pets.size() <= 1) {
            return distances;
        }

        Collections.sort(pets);

        for (int i = 0; i < pets.size(); i++) {
            PET pet = pets.elementAt(i);
            int minDistance = Integer.MAX_VALUE;
            // search descending
            for (int j = i - 1; j >= 0; j--) {
                if (HIT2.calculateDistance(pet.getHead(), pets.elementAt(j).getHead()) > minDistance) {
                    break;
                }
                int distance1 = HIT2.calculateDistance(pet.getHead(), pets.elementAt(j).getHead());
                int distance2 = HIT2.calculateDistance(pet.getTail(), pets.elementAt(j).getTail());
                if (distance1 > Integer.MAX_VALUE / 3) {
                    distance1 = Integer.MAX_VALUE / 3;
                }
                if (distance2 > Integer.MAX_VALUE / 3) {
                    distance2 = Integer.MAX_VALUE / 3;
                }

                int distance = distance1 + distance2;
                if (minDistance > distance) {
                    minDistance = distance;
                }

            }
            // search ascending
            for (int j = i + 1; j < pets.size(); j++) {
                if (HIT2.calculateDistance(pet.getHead(), pets.elementAt(j).getHead()) > minDistance) {
                    break;
                }
                int distance1 = HIT2.calculateDistance(pet.getHead(), pets.elementAt(j).getHead());
                int distance2 = HIT2.calculateDistance(pet.getTail(), pets.elementAt(j).getTail());
                if (distance1 > Integer.MAX_VALUE / 3) {
                    distance1 = Integer.MAX_VALUE / 3;
                }
                if (distance2 > Integer.MAX_VALUE / 3) {
                    distance2 = Integer.MAX_VALUE / 3;
                }

                int distance = distance1 + distance2;
                if (minDistance > distance) {
                    minDistance = distance;
                }
            }

            distances.add(new Integer(minDistance));
        }
        return distances;
    }

    public static Vector<PET> removeWeakPets(Vector<PET> pets, double weightCutoff) {
        if (pets.size() <= 1) {
            return pets;
        }

        Vector<PET> newChiapets = new Vector<PET>(); // Interaction-PET
        for (int i = 0; i < pets.size(); i++) {
            if (pets.elementAt(i).getWeight() >= weightCutoff) {
                newChiapets.add(pets.elementAt(i));
            }
        }
        return newChiapets;
    }

    /**
     * @return the index
     */
    public int getIndex() {
        return index;
    }

    /**
     * @param index the index to set
     */
    public void setIndex(int index) {
        this.index = index;
    }
}
