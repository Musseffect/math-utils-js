

export default abstract class IConstMatrix {
    abstract get(row: number, col: number): number;
    abstract numRows(): number;
    abstract numCols(): number;
    isSquare(): boolean {
        return this.numCols() == this.numRows();
    }
    height(): number {
        return this.numRows();
    }
    width(): number {
        return this.numCols();
    }
}