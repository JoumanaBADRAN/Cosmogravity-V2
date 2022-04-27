/**
 * @class Graphic
 *
 * attributes :
 * @param id_graph Identifier
 * @param simulation association between graphic and simulation
 */
export class Graphic {
    //-------------------------constructor-----------------------
    constructor(id_graph, simulation) {
        this._id_graph = id_graph;
        this._simulation = simulation;
    }
    //--------------------------Accessors------------------------
    // id_graph
    get id_graph() {
        return this._id_graph;
    }
    // simulation
    get simulation() {
        return this._simulation;
    }
}
