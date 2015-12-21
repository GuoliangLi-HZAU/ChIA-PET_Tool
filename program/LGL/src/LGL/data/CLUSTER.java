/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

/**
 *
 * The head and tail in a cluster are sorted by their genomic locations.
 *
 * @author ligl
 */
import LGL.util.SeqUtil;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

public class CLUSTER implements Comparable {

    private REGION head;
    private REGION tail;
    private double weight = 0.0;
    private String annotation;
    private boolean bSearched;
    private String petIndexes = "---";

    public CLUSTER() {
        head = null;
        tail = null;
        this.setAnnotation("---");
        this.setbSearched(false);
    }

    public CLUSTER(PET pet, int extensionLength) {
        if (SeqUtil.isForwardStrand(pet.getHead().getStrand()) == true) {
            head = new REGION(pet.getHead().getChrom(), pet.getHead().getLoci(), pet.getHead().getLoci() + extensionLength);
        } else {
            head = new REGION(pet.getHead().getChrom(), pet.getHead().getLoci() - extensionLength, pet.getHead().getLoci());
        }

        if (SeqUtil.isForwardStrand(pet.getTail().getStrand()) == true) {
            tail = new REGION(pet.getTail().getChrom(), pet.getTail().getLoci() - extensionLength, pet.getTail().getLoci());
        } else {
            tail = new REGION(pet.getTail().getChrom(), pet.getTail().getLoci(), pet.getTail().getLoci() + extensionLength);
        }
        this.weight = 1.0;
        this.setAnnotation("---");
        this.setbSearched(false);
        sortHeadTail();
    }

    public CLUSTER(PET pet, int extensionLength, int sortingLabel, Hashtable<String, Integer> chromLength) {
        head = null;
        tail = null;
        this.weight = 1.0;
        this.setAnnotation("---");
        this.setbSearched(false);

        if (chromLength.containsKey(pet.getHead().getChrom())) {
            if (SeqUtil.isForwardStrand(pet.getHead().getStrand()) == true) {
                if (pet.getHead().getLoci() - extensionLength > 0) {
                    head = new REGION(pet.getHead().getChrom(), pet.getHead().getLoci() - extensionLength, pet.getHead().getLoci());
                } else {
                    head = new REGION(pet.getHead().getChrom(), 1, pet.getHead().getLoci());
                }
            } else {
                if (pet.getHead().getLoci() + extensionLength < chromLength.get(pet.getHead().getChrom())) {
                    head = new REGION(pet.getHead().getChrom(), pet.getHead().getLoci(), pet.getHead().getLoci() + extensionLength);
                } else {
                    head = new REGION(pet.getHead().getChrom(), pet.getHead().getLoci(), chromLength.get(pet.getHead().getChrom()));
                }
            }
        }
        if (chromLength.containsKey(pet.getTail().getChrom())) {
            if (SeqUtil.isForwardStrand(pet.getTail().getStrand()) == true) {
                if (pet.getTail().getLoci() - extensionLength > 0) {
                    tail = new REGION(pet.getTail().getChrom(), pet.getTail().getLoci() - extensionLength, pet.getTail().getLoci());
                } else {
                    tail = new REGION(pet.getTail().getChrom(), 1, pet.getTail().getLoci());
                }
            } else {

                if (pet.getTail().getLoci() + extensionLength < chromLength.get(pet.getTail().getChrom())) {
                    tail = new REGION(pet.getTail().getChrom(), pet.getTail().getLoci(), pet.getTail().getLoci() + extensionLength);
                } else {
                    tail = new REGION(pet.getTail().getChrom(), pet.getTail().getLoci(), chromLength.get(pet.getTail().getChrom()));
                }
            }
        }

        if (sortingLabel == 1) {
            sortHeadTail();
        } else if (sortingLabel == -1) {
            sortHeadTail_descending();
        }
    }

    public CLUSTER(REGION head, REGION tail) {
        this.head = head;
        this.tail = tail;
        this.weight = 1.0;
        this.setAnnotation("---");
        this.setbSearched(false);
        sortHeadTail();
    }

    public CLUSTER(REGION head, REGION tail, double weight) {
        this.head = head;
        this.tail = tail;
        this.weight = weight;
        this.setAnnotation("---");
        this.setbSearched(false);
        sortHeadTail();
    }

    public CLUSTER(REGION head, REGION tail, double weight, int sortingLabel) {
        this.head = head;
        this.tail = tail;
        this.weight = weight;
        this.setAnnotation("---");
        this.setbSearched(false);

        if (sortingLabel == 1) {
            sortHeadTail();
        } else if (sortingLabel == -1) {
            sortHeadTail_descending();
        }
    }

    public CLUSTER(REGION head, REGION tail, double weight, String annotation) {
        this.head = head;
        this.tail = tail;
        this.weight = weight;
        this.setAnnotation(annotation);
        this.setbSearched(false);
        sortHeadTail();
    }

    private void sortHeadTail() {
        if (this.getHead().compareTo(this.getTail()) > 0) {
            REGION tempRegion = this.getHead();
            this.setHead(this.getTail());
            this.setTail(tempRegion);
        }
    }

    private void sortHeadTail_descending() {
        if (this.getHead().compareTo(this.getTail()) < 0) {
            REGION tempRegion = this.getHead();
            this.setHead(this.getTail());
            this.setTail(tempRegion);
        }
    }

    @Override
    public String toString() {
        return (this.getHead().toString() + "\t" + this.getTail().toString() + "\t" + this.getWeight() + "\t" + this.getPetIndexes());
    }

    public String toString(String sep) {
        return (this.getHead().toString(sep) + sep + this.getTail().toString(sep) + sep + this.getWeight() + sep + this.getPetIndexes());
    }

    public String toString(String sep, int mode) {
        return (this.getHead().toString(mode) + sep + this.getTail().toString(mode) + sep + this.getWeight() + sep + this.getPetIndexes());
    }

    public String toString2() {
        int sameChrom = 0;
        if (this.getHead().getChrom().compareTo(this.getTail().getChrom()) == 0) {
            sameChrom = 1;
        }
        return (this.getHead().toString() + "\t" + this.getTail().toString() + "\t" + sameChrom + "\t" + REGION.calculateCenterDistance(this.getHead(), this.getTail()) + "\t" + this.getWeight() + "\t" + this.getPetIndexes());
    }

    public void combine(CLUSTER anotherCluster) throws Exception {
        this.getHead().combine(anotherCluster.getHead());
        this.getTail().combine(anotherCluster.getTail());
        this.setWeight(this.getWeight() + anotherCluster.getWeight());
    }

    public int compareTo(Object anotherCluster) throws ClassCastException {
        if (!(anotherCluster instanceof CLUSTER)) {
            throw new ClassCastException("A CLUSTER object expected.");
        }
        int result = this.getHead().compareTo(((CLUSTER) anotherCluster).getHead());
        if (result == 0) {
            result = this.getTail().compareTo(((CLUSTER) anotherCluster).getTail());
        }
        return result;
    }

    /**
     * @return the head
     */
    public REGION getHead() {
        return head;
    }

    /**
     * @param head the head to set
     */
    public void setHead(REGION head) {
        this.head = head;
    }

    /**
     * @return the tail
     */
    public REGION getTail() {
        return tail;
    }

    /**
     * @param tail the tail to set
     */
    public void setTail(REGION tail) {
        this.tail = tail;
    }

    public boolean overlap(CLUSTER anotherCluster) {
        boolean overlapped = false;
        if (((this.getHead().overlap(anotherCluster.getHead()) == true) && (this.getTail().overlap(anotherCluster.getTail()) == true)) || ((this.getHead().overlap(anotherCluster.getTail()) == true) && (this.getTail().overlap(anotherCluster.getHead()) == true))) {
            overlapped = true;
        }
        return overlapped;
    }

    public boolean overlap_orderSpecific(CLUSTER anotherCluster) {
        boolean overlapped = false;
        if ((this.getHead().overlap(anotherCluster.getHead()) == true) && (this.getTail().overlap(anotherCluster.getTail()) == true)) {
            overlapped = true;
        }
        return overlapped;
    }

    public boolean overlap(CLUSTER2 cluster2) {
        boolean overlapped = false;
        if (((this.getHead().overlap(cluster2.getHead()) == true) && (this.getTail().overlap(cluster2.getTail()) == true)) || ((this.getHead().overlap(cluster2.getTail()) == true) && (this.getTail().overlap(cluster2.getHead()) == true))) {
            overlapped = true;
        }
        return overlapped;
    }

    public int calculateClusterAnchorBoundaryDistance(CLUSTER anotherCluster) {
        if (this.overlap(anotherCluster)) {
            return 0; // if the clusters are overlapped at both anchors, the cluster anchor boundary distance is 0
        }

        int distance = Integer.MAX_VALUE;
        int distStart1_Start2 = REGION.calculateBoundaryDistance(this.getHead(), anotherCluster.getHead());
        int distStart1_End2 = REGION.calculateBoundaryDistance(this.getHead(), anotherCluster.getTail());
        int distEnd1_Start2 = REGION.calculateBoundaryDistance(this.getTail(), anotherCluster.getHead());
        int distEnd1_End2 = REGION.calculateBoundaryDistance(this.getTail(), anotherCluster.getTail());

        if ((distStart1_Start2 < Integer.MAX_VALUE) && (distEnd1_End2 < Integer.MAX_VALUE)) {
            if (distance > distStart1_Start2 + distEnd1_End2) {
                distance = distStart1_Start2 + distEnd1_End2;
            }
        }
        if ((distStart1_End2 < Integer.MAX_VALUE) && (distEnd1_Start2 < Integer.MAX_VALUE)) {
            if (distance > distStart1_End2 + distEnd1_Start2) {
                distance = distStart1_End2 + distEnd1_Start2;
            }
        }

        return distance;
    }

    public static Vector<REGION3> cluster2Region3(CLUSTER cluster) {
        Vector<REGION3> region3s = new Vector<REGION3>();
        REGION3 region3Head = new REGION3(cluster.getHead(), (int) cluster.getWeight());
        REGION3 region3Tail = new REGION3(cluster.getTail(), (int) cluster.getWeight());
        region3s.add(region3Head);
        region3s.add(region3Tail);
        return region3s;
    }

    // convert clusters to Region3 format
    public static Vector<REGION3> cluster2Region3(Vector<CLUSTER> clusters) {
        Vector<REGION3> region3s = new Vector<REGION3>();
        for (int iCluster = 0; iCluster < clusters.size(); iCluster++) {
            REGION3 region3Head = new REGION3(clusters.elementAt(iCluster).getHead(), (int) clusters.elementAt(iCluster).getWeight());
            REGION3 region3Tail = new REGION3(clusters.elementAt(iCluster).getTail(), (int) clusters.elementAt(iCluster).getWeight());
            region3s.add(region3Head);
            region3s.add(region3Tail);
        }
        return region3s;
    }

    public static Vector<CLUSTER> load(String ClusterFile) throws IOException {
        Vector<CLUSTER> clusters = new Vector<CLUSTER>();
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(ClusterFile))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String fields[] = line.split("\t");
            REGION head = new REGION(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]));
            REGION tail = new REGION(fields[3], Integer.parseInt(fields[4]), Integer.parseInt(fields[5]));

            double weight = 1.0;
            if (fields.length >= 7) {
                weight = Double.parseDouble(fields[6]);
            }
            String petIndexes = "---";
            if (fields.length >= 8) {
                petIndexes = fields[7];
            }

            CLUSTER cluster = new CLUSTER(head, tail, weight);
            cluster.setPetIndexes(petIndexes);
            clusters.add(cluster);
        }
        fileIn.close();

        //Collections.sort(clusters);
        return clusters;
    }

    // different from load()
    // "OrderSpecific" means that the head and tail regions are kept in the same order
    public static Vector<CLUSTER> load_orderSpecific(String ClusterFile) throws IOException {
        int sortingLabel = 0; // no sorting on the head and tail regions
        Vector<CLUSTER> clusters = new Vector<CLUSTER>();
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(ClusterFile))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String fields[] = line.split("\t");
            REGION head = new REGION(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]));
            REGION tail = new REGION(fields[3], Integer.parseInt(fields[4]), Integer.parseInt(fields[5]));

            double weight = 1.0;
            if (fields.length >= 7) {
                weight = Double.parseDouble(fields[6]);
            }
            String petIndexes = "---";
            if (fields.length >= 8) {
                petIndexes = fields[7];
            }

            CLUSTER cluster = new CLUSTER(head, tail, weight, sortingLabel);
            cluster.setPetIndexes(petIndexes);
            clusters.add(cluster);
        }
        fileIn.close();

        //Collections.sort(clusters);
        return clusters;
    }

    public static int getMaxRegionSpan(Vector<CLUSTER> clusters) {
        int maxRegionSpan = 0;
        for (int i = 0; i < clusters.size(); i++) {
            int span = clusters.elementAt(i).getHead().getSpan();
            if (maxRegionSpan < span) {
                maxRegionSpan = span;
            }
            span = clusters.elementAt(i).getTail().getSpan();
            if (maxRegionSpan < span) {
                maxRegionSpan = span;
            }
        }
        return maxRegionSpan;
    }

    public static void saveClusters(String ClusterFile, Vector<CLUSTER> clusters) throws IOException {
        PrintWriter fileOut = new PrintWriter(new BufferedWriter(new FileWriter(ClusterFile)));

        for (int i = 0; i < clusters.size(); i++) {
            fileOut.println(clusters.elementAt(i).toString());
        }

        fileOut.close();
    }

    public static Hashtable<String, Vector<CLUSTER>> cluster2Hash(Vector<CLUSTER> clusters) {
        // separate the clusters by head chromosomes
        Hashtable<String, Vector<CLUSTER>> cluster2Hash = new Hashtable<String, Vector<CLUSTER>>();
        for (int i = 0; i < clusters.size(); i++) {
            Vector<CLUSTER> clustersTemp = cluster2Hash.get(clusters.elementAt(i).getHead().getChrom());
            if (clustersTemp == null) {
                clustersTemp = new Vector<CLUSTER>();
                cluster2Hash.put(clusters.elementAt(i).getHead().getChrom(), clustersTemp);
            }
            clustersTemp.add(clusters.elementAt(i));
        }
        // sort the regions in each chromosome
        Set<String> chromNames = cluster2Hash.keySet();
        if (chromNames != null) {
            Iterator<String> itr = chromNames.iterator();
            while (itr.hasNext()) {
                String chromName = itr.next();
                Vector<CLUSTER> clustersTemp = (Vector<CLUSTER>) cluster2Hash.get(chromName);
                Collections.sort(clustersTemp);
                cluster2Hash.put(chromName, clustersTemp);
                //System.out.println("chromName" + chromName + ":\t" +clustersTemp.size() );
            }
        }

        return cluster2Hash;
    }

    public static boolean clusterOverlap(CLUSTER cluster, Hashtable<String, Vector<CLUSTER>> cluster2Hash, int maxRegionSpan) {
        boolean overlapped = false;

        Vector clustersTemp = cluster2Hash.get(cluster.getHead().getChrom());
        if (clustersTemp != null) {
            int index = Collections.binarySearch(clustersTemp, cluster);
            if (index < 0) {
                index = -index - 1;
            }
            if (index >= clustersTemp.size()) {
                index = clustersTemp.size() - 1;
            }

            if (cluster.overlap(((Vector<CLUSTER>) clustersTemp).elementAt(index)) == true) {
                overlapped = true;
                return overlapped;
            }

            for (int iCluster = index - 1; iCluster >= 0; iCluster--) {
                // exhaustive search, in case the regions are overlapped
                // if the regions are not overlapped, binary search can be faster
                if (cluster.overlap(((Vector<CLUSTER>) clustersTemp).elementAt(iCluster)) == true) {
                    overlapped = true;
                    return overlapped;
                } else if (cluster.getHead().getStart() - ((Vector<CLUSTER>) clustersTemp).elementAt(iCluster).getHead().getStart() > maxRegionSpan) {
                    break;
                }
//                } else if (cluster.getHead().overlap(((Vector<CLUSTER>) clustersTemp).elementAt(iCluster).getHead()) == false) {
//                    break;
//                }
            }
            for (int iCluster = index + 1; iCluster < clustersTemp.size(); iCluster++) {
                // exhaustive search, in case the regions are overlapped
                // if the regions are not overlapped, binary search can be faster
                if (cluster.overlap(((Vector<CLUSTER>) clustersTemp).elementAt(iCluster)) == true) {
                    overlapped = true;
                    return overlapped;
                } else if (cluster.getHead().overlap(((Vector<CLUSTER>) clustersTemp).elementAt(iCluster).getHead()) == false) {
                    break;
                }
            }
        }
        return overlapped;
    }

    public static void clusterOverlap(CLUSTER cluster, Hashtable<String, Vector<CLUSTER>> cluster2Hash, int maxRegionSpan, PrintWriter fileOut) throws IOException {
        fileOut.print(cluster.toString());
        boolean firstOverlap = true;

        Vector clustersTemp = cluster2Hash.get(cluster.getHead().getChrom());
        //System.out.println("chromName 1: " + cluster.getHead().getChrom());
        if (clustersTemp != null) {
            int index = Collections.binarySearch(clustersTemp, cluster);
            //System.out.println("index = " + index);
            if (index < 0) {
                index = -index - 1;
            }
            if (index >= clustersTemp.size()) {
                index = clustersTemp.size() - 1;
            }
            //System.out.println(cluster.toString() + " <==> " + ((Vector<CLUSTER>) clustersTemp).elementAt(index).toString()) ;

            if (cluster.overlap(((Vector<CLUSTER>) clustersTemp).elementAt(index)) == true) {
                if (firstOverlap == true) {
                    fileOut.print("\tYes");
                    firstOverlap = false;
                }
                fileOut.print("\t" + ((Vector<CLUSTER>) clustersTemp).elementAt(index).toString());
            }

            for (int iCluster = index - 1; iCluster >= 0; iCluster--) {
                // exhaustive search, in case the regions are overlapped
                // if the regions are not overlapped, binary search can be faster
                if (cluster.overlap(((Vector<CLUSTER>) clustersTemp).elementAt(iCluster)) == true) {
                    if (firstOverlap == true) {
                        fileOut.print("\tYes");
                        firstOverlap = false;
                    }
                    fileOut.print("\t" + ((Vector<CLUSTER>) clustersTemp).elementAt(iCluster).toString());
                } else if (cluster.getHead().getStart() - ((Vector<CLUSTER>) clustersTemp).elementAt(iCluster).getHead().getStart() > maxRegionSpan) {
                    break;
                }
            }
            for (int iCluster = index + 1; iCluster < clustersTemp.size(); iCluster++) {
                // exhaustive search, in case the regions are overlapped
                // if the regions are not overlapped, binary search can be faster
                if (cluster.overlap(((Vector<CLUSTER>) clustersTemp).elementAt(iCluster)) == true) {
                    if (firstOverlap == true) {
                        fileOut.print("\tYes");
                        firstOverlap = false;
                    }
                    fileOut.print("\t" + ((Vector<CLUSTER>) clustersTemp).elementAt(iCluster).toString());
                } else if (cluster.getHead().overlap(((Vector<CLUSTER>) clustersTemp).elementAt(iCluster).getHead()) == false) {
                    break;
                }
            }
        }
        if (firstOverlap == true) {
            fileOut.println("\tNo");
        } else {
            fileOut.println();
        }
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

    /**
     * @return the annotation
     */
    public String getAnnotation() {
        return annotation;
    }

    /**
     * @param annotation the annotation to set
     */
    public void setAnnotation(String annotation) {
        this.annotation = annotation;
    }

    /**
     * @return the bSearched
     */
    public boolean isbSearched() {
        return bSearched;
    }

    /**
     * @param bSearched the bSearched to set
     */
    public void setbSearched(boolean bSearched) {
        this.bSearched = bSearched;
    }

    /**
     * @return the petIndexes
     */
    public String getPetIndexes() {
        return petIndexes;
    }

    /**
     * @param petIndexes the petIndexes to set
     */
    public void setPetIndexes(String petIndexes) {
        this.petIndexes = petIndexes;
    }
}
