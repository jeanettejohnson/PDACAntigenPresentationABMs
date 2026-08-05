from omeify.inputs import AkoyaMIFQptiff, AkoyaHEQptiff

input_processor = AkoyaMIFQptiff(input_file_path, series=series_number)
input_processor.rename_channels = rename_channels_dict

output_info = input_processor.convert(output_file_path, display_uuid=True)