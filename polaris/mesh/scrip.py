from mpas_tools.logging import check_call
from pyremap import MpasCellMeshDescriptor


def write_mpas_scrip_file(mesh_filename, scrip_filename, mesh_name, logger):
    """
    Write a SCRIP file describing the cells of an MPAS mesh.

    Parameters
    ----------
    mesh_filename : str
        The MPAS mesh file to describe

    scrip_filename : str
        The SCRIP file to write, in NETCDF3_64BIT_DATA format

    mesh_name : str
        The name to give the mesh in the SCRIP file

    logger : logging.Logger
        A logger for the output from ``ncks``
    """
    netcdf4_filename = scrip_filename.replace('.nc', '.netcdf4.nc')

    descriptor = MpasCellMeshDescriptor(
        filename=mesh_filename,
        mesh_name=mesh_name,
    )
    descriptor.to_scrip(netcdf4_filename)

    # writing directly in NETCDF3_64BIT_DATA proved prohibitively slow
    # so we will use ncks to convert
    args = [
        'ncks',
        '-O',
        '-5',
        netcdf4_filename,
        scrip_filename,
    ]
    check_call(args, logger)
